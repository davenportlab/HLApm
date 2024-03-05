#!/usr/bin/env python2

import collections
import gzip
import operator
import optparse
import os
import re
import sys
import time

import pybam


# Define a wrapper around quicksect.IntervalTree or intervaltree.IntervalTree, depending
# on which one is available, so that the rest of the script can work with both.
# quicksect is much faster than intervaltree, but harder to install.

# Throughout the script, start and end coordinates are all 1-based and both included
try:
    import quicksect

    class myintervaltree(quicksect.IntervalTree):

        def __init__(self):
            super(myintervaltree, self).__init__()
            self.interval_to_data = {}
            self.has_overlap = self.search

        # quicksect cannot associate intervals to data, so we need to it ourselves
        def add_data(self, start, end, data):
            interval = quicksect.Interval(start, end)
            self.insert( interval )
            self.interval_to_data[ interval ] = data

        def data_overlapping(self, start, end):
            return [self.interval_to_data[i] for i in self.search(start, end)]

except ImportError:
    import intervaltree

    class myintervaltree(intervaltree.IntervalTree):

        def __init__(self):
            super(myintervaltree, self).__init__()
            # Deal with versions 2.* and 3.*
            if hasattr(intervaltree.IntervalTree, "search"):
                self.search_func = self.search
            else:
                self.search_func = self.overlap
            self.has_overlap = self.search_func

        # In intervaltree.IntervalTree the intervals include the lower bound but not the upper bound.
        def add_data(self, start, end, data):
            self.addi(start, end+1, data)

        def has_overlap(self, start, end):
            return self.search_func(start, end+1)

        def data_overlapping(self, start, end):
            return [i.data for i in self.search_func(start, end+1)]


# Simple structure for a read alignment
ReadAlignment = collections.namedtuple('ReadAlignment', ['chrom_name', 'start_pos', 'cigar_list', 'NH_value', 'NM_value'])
# Simple structure for a generic, named, data
NamedEntry = collections.namedtuple('NamedEntry', ['name', 'data'])


parser = optparse.OptionParser(usage = "usage: %prog [options] gtf_file bam_file_1 bam_file_2 ...")
parser.add_option("-r", "--rename_gene", nargs = 2, action = 'append', dest = 'gene_renames', metavar = 'NAME_TO_REPLACE NEW_NAME',
        help = 'Replace some erroneous gene names')
parser.add_option("-g", "--gene_types", dest = 'gene_types', default = 'protein_coding',
        help = 'Comma-separated list of gene biotypes to use [default: %default]. Use an empty string for no filtering')
parser.add_option("-t", "--transcript_types", dest = 'transcript_types', default = 'protein_coding',
        help = 'Comma-separated list of transcript biotypes to use for the exon-overlap filtering [default: %default]. Use an empty string for no filtering')
parser.add_option("-n", "--no_gtf_filter", action = 'store_true',
        help = 'Do not use a GTF file to filter the reads. The command line arguments are then expected to all be BAM files.')
parser.add_option("-s", "--strict", action = 'store_true', default = False,
        help = 'Use to select the strict mode, which discards paired reads when one read alf is not in an exon')
parser.add_option("-l", "--log", nargs = 1, dest = 'log_file', action = 'store',
        help = "Do not include paired reads when one read is not in an exon")

(options, args) = parser.parse_args()

if options.no_gtf_filter:
    parser.usage = parser.usage.replace('gtf_file ', '')
    gtf_file = None
else:
    if len(args) == 0:
        parser.error("No GTF/BAM files given")
    gtf_file = args.pop(0)

if len(args) == 0:
    parser.error("No BAM files given")

bam_files = args
print >> sys.stderr, "BAM+GTF merger for AltHapAlign called with these options"
print >> sys.stderr, "\tgtf_file:", gtf_file
print >> sys.stderr, "\tbam_files:", bam_files
print >> sys.stderr, "\toptions:", options
print >> sys.stderr

class StopWatch:
    def __init__(self):
        self.ref_time = time.time()

    def reset(self):
        old_time = self.ref_time
        self.ref_time = time.time()
        return self.ref_time-old_time

main_stopwatch = StopWatch()

def load_genes(gtf_file):
    if not options.no_gtf_filter:
        print >> sys.stderr, "Loading the GTF file ... ",
    gene_type_filter = re.compile('gene_type "?(%s)"?;' % options.gene_types.replace(",", "|")) if options.gene_types else None
    transcript_type_filter = re.compile('transcript_type "?(%s)"?;' % options.transcript_types.replace(",", "|")) if options.transcript_types else None
    gene_renames = dict(options.gene_renames) if options.gene_renames else {}

    gene_names = collections.defaultdict(myintervaltree)
    n_genes = set()
    n_exons = 0
    if options.no_gtf_filter:
        gtf_file = os.devnull
    with gzip.GzipFile(gtf_file, 'r') if gtf_file.endswith('.gz') else open(gtf_file, 'r') as f:
        for line in f:
            t = line.strip().split("\t")
            if gene_type_filter and not gene_type_filter.search(t[8]):
                continue
            if transcript_type_filter and not transcript_type_filter.search(t[8]):
                continue
            if t[2] == "exon":
                if 'gene_name' not in t[8]:
                    print >> sys.stderr, "Failed\nERROR: no 'gene_name' attribute for this exon:", line,
                    sys.exit(1)
                i1 = t[8].find('gene_name')
                i2 = t[8].find(';', i1)
                gn = t[8][i1+9:i2].strip().strip('"')
                name = gene_renames.get(gn, gn)
                start = int(t[3])
                end = int(t[4])
                gene_names[t[0]].add_data(start, end, name)
                n_genes.add(name)
                n_exons += 1
    if not options.no_gtf_filter:
        print >> sys.stderr, "Done (%.2f seconds): %d genes and %d exons" % (main_stopwatch.reset(), len(n_genes), n_exons)
    return gene_names

gene_names = load_genes(gtf_file)

# Input: BAM read with tag list, and tag name
# Output: the value of the tag
# Exception: exit if the tag is not present, or in multiple copies
def get_tag_value(read, tag):
    values = [t[2] for t in read[-1] if t[0].upper() == tag]
    if not values:
        print >> sys.stderr, "ERROR: Could not find any %s tag in" % tag, read
        sys.exit(1)
    if len(values) > 1:
        print >> sys.stderr, "ERROR: More than 1 %s tag in" % tag, read
        sys.exit(1)
    return values[0]


# Input: CIGAR alignment already parsed to a list of (number, character)
# Output: integer
# Description: Computes the length of the alignment on the genome
def mapping_length(cigar_list):
    s = 0
    for (n, c) in cigar_list:
        if c == 'I':
            s -= n
        else:
            s += n
    return s


# Holds all the BAM processing functions (generators) and binds them together
class BAMProcessor:

    def __init__(self, bam_file):
        self.bam_parser = pybam.read(bam_file, ['sam_qname', 'sam_rname', 'sam_pos1', 'sam_cigar_list', 'sam_tags_list'])
        self.n_bam_aligns = 0
        self.n_singletons = 0
        self.n_paired_alignments = 0
        self.n_multiple_hits = 0
        self.n_non_exonic_bam_aligns = 0
        self.n_different_genes_pair = 0
        self.n_partially_exonic = 0
        self.n_fully_exonic = 0
        self.n_gene_match = 0


    # Input: BAM iterator with tag list
    # Output: BAM iterator with selected tag values (NH and NM)
    # Description: Extract the NH and NM values from the BAM alignments
    def extract_tags(self):
        for read in self.bam_parser:
            self.n_bam_aligns += 1
            NH_value = get_tag_value(read, 'NH')
            NM_value = get_tag_value(read, 'NM')
            yield NamedEntry(read[0], ReadAlignment(read[1], read[2], read[3], NH_value, NM_value))


    # Input: iterator (read_name, ...data...)
    # Output: iterator (read_name, [(...data1...), (...data2...)])
    # Description: Group consecutive pairs of reads that have the same name.
    #              This assumes that the input BAM file is sorted.
    def group_read_alignments(self):
        last_reads = []
        last_read_name = None
        for read in self.extract_tags():
            this_read_name = read.name
            if last_read_name:
                if last_read_name != this_read_name:
                    yield NamedEntry(last_read_name, last_reads)
                    last_reads = []
                    last_read_name = this_read_name
            else:
                last_read_name = this_read_name
            last_reads.append(read.data)
        if last_read_name:
            yield NamedEntry(last_read_name, last_reads)


    # Input: iterator over paired alignments
    # Output: iterator over a subset of the input
    # Description: Discards the alignments that are not uniquely mapped pairs,
    #              i.e. are singletons or have multiple hits / don't have NH:i:1.
    def select_paired_alignments(self):
        for aligned_read in self.group_read_alignments():
            n_alignments = len(aligned_read.data)
            NH_values = [alignment.NH_value for alignment in aligned_read.data]
            if n_alignments == 1:
                self.n_singletons += 1
            elif n_alignments == 2:
                # Check that NH = 1 for both
                if NH_values[0] != 1:
                    assert NH_values[1] != 1
                    # Could be classified as singletons too
                    self.n_multiple_hits += 2
                else:
                    assert NH_values[1] == 1
                    self.n_paired_alignments += 1
                    yield aligned_read
            else:
                for NH_value in NH_values:
                    assert NH_value != 1
                self.n_multiple_hits += n_alignments


    # Input: iterator over paired alignments
    # Output: iterator that has gene names too, but on a subset of the alignments
    # Description: Discards the reads that do not overlap any exons
    def only_exonic_mappings(self):
        for aligned_read in self.select_paired_alignments():
            ok = 0
            genes_seen = set()
            for aln in aligned_read.data:
                if aln.chrom_name in gene_names:
                    end_pos = aln.start_pos + mapping_length(aln.cigar_list) - 1
                    names_here = gene_names[aln.chrom_name].data_overlapping(aln.start_pos, end_pos)
                    if names_here:
                        genes_seen.update(names_here)
                        ok += 1
            if ok == 0:
                self.n_non_exonic_bam_aligns += 1
            elif len(genes_seen) > 1:
                self.n_different_genes_pair += 1
            else:
                if ok == 1:
                    self.n_partially_exonic += 1
                    if options.strict:
                        continue
                else:
                    self.n_fully_exonic += 1
                yield NamedEntry(aligned_read.name, NamedEntry(genes_seen.pop(), aligned_read.data))

class BAMMergedProcessor:
    def __init__(self, bam_processors, generator_name):
        self.bam_processors = bam_processors
        self.bam_iterators = map(operator.methodcaller(generator_name), bam_processors)
        self.n_groups = 0
        self.n_unique_groups = 0
        self.n_best_groups = 0
        self.n_ambiguous_groups = 0

    def stat_sum(self, attr_name):
        return sum(map(operator.attrgetter(attr_name), self.bam_processors))


    # Input: list of iterators (read_name, data)
    # Output: iterator (read_name, [data1 | None, data2 | None, ...])
    # Description: For each read name, groups the reads from each BAM file,
    #              assuming that the files are sorted by read name
    # Remarks: This is similar to a k-way merge algorithm. Although the merge
    #          itself can be achieved in O(log(k)) in theory, my implementations
    #          did not provide any improvements, especially because the output
    #          of the function is O(k). So in the end I'm merging the streams in O(k)
    def merged_iterators(self):
        EMPTY_ITERATOR = NamedEntry(None, None)
        current_reads = []
        active_bam_parsers = 0
        for parser in self.bam_iterators:
            try:
                current_reads.append( parser.next() )
                active_bam_parsers += 1
            except StopIteration:
                current_reads.append( EMPTY_ITERATOR )
        last_read_name = None
        while active_bam_parsers:
            read_names = [cr.name for cr in current_reads if cr is not EMPTY_ITERATOR]
            next_read_name = min(read_names)
            if last_read_name and (next_read_name < last_read_name):
                not_sorted = []
                for (i, cr) in enumerate(current_reads):
                    if cr.name == next_read_name:
                        not_sorted.append(bam_files[i])
                raise AssertionError("The BAM files %s are not sorted: read '%s' is after '%s'" % (", ".join(not_sorted), next_read_name, last_read_name))
            last_read_name = next_read_name
            self.n_groups += 1
            yield NamedEntry(next_read_name, [cr.data if cr.name == next_read_name else None for cr in current_reads])
            for (i, cr) in enumerate(current_reads):
                if cr.name == next_read_name:
                    try:
                        current_reads[i] = self.bam_iterators[i].next()
                    except StopIteration:
                        current_reads[i] = EMPTY_ITERATOR
                        active_bam_parsers -= 1


    # Input: iterator (read_name, [data1 | None, data2 | None, ...])
    # Output: iterator (read_name, [data1 | None, data2 | None, ...])
    # Description: For each group of reads, finds the gene names on the genome
    #              (using the alignment) and only keeps the groups that map to
    #              a single gene name.
    def select_same_gene(self):
        for read_group in self.merged_iterators():
            genes_seen = {}
            for pair in read_group.data:
                if pair is not None:
                    gene_name = pair.name
                    NM_score = sum(alignment.NM_value for alignment in pair.data)
                    if gene_name in genes_seen:
                        if NM_score < genes_seen[gene_name]:
                            genes_seen[gene_name] = NM_score
                    else:
                        genes_seen[gene_name] = NM_score
            if len(genes_seen) == 1:
                status = 'unique'
                best_genes = genes_seen.keys()
                self.n_unique_groups += 1
            else:
                best_NM = min(genes_seen.values())
                best_genes = [gene_name for (gene_name, NM_score) in genes_seen.items() if NM_score == best_NM]
                if len(best_genes) == 1:
                    self.n_best_groups += 1
                    status = 'best'
                else:
                    self.n_ambiguous_groups += 1
                    status = 'ambiguous'
            for gene_name in best_genes:
                for (bp, pair) in zip(self.bam_processors, read_group.data):
                    if pair and pair.name == gene_name:
                        bp.n_gene_match += 1
                yield NamedEntry(read_group.name, (status, gene_name, [pair.data if pair and pair.name == gene_name else None for pair in read_group.data]))


def toString(entry):
    line = [entry.name]
    if isinstance(entry.data, tuple):
        # gene name
        line.extend(entry.data[:-1])
        alns = entry.data[-1]
    else:
        alns = entry.data
    for paired_data in alns:
        if paired_data is None:
            line.append('NA')
        else:
            line.append(str(paired_data[0].NM_value + paired_data[1].NM_value))
    return "\t".join(line)


print >> sys.stderr, "Opening the BAM files ...",
bam_processors = map(BAMProcessor, bam_files)
print >> sys.stderr, " Done (%.2f seconds)" % main_stopwatch.reset()

last_n_bam_aligns = 0
print >> sys.stderr, "Reading the BAM files ..."
headers = ["read_name"]
if options.no_gtf_filter:
    merged_processor = BAMMergedProcessor(bam_processors, 'select_paired_alignments')
    it = merged_processor.merged_iterators()
else:
    headers = headers + ["gene_name_confidence", "gene_name"]
    merged_processor = BAMMergedProcessor(bam_processors, 'only_exonic_mappings')
    it = merged_processor.select_same_gene()
headers = headers + bam_files
print "\t".join(headers)
last_n_groups = 0
n_lines = 0
bam_reading_stopwatch = StopWatch()
for entry in it:
    print toString(entry)
    n_lines += 1
    if merged_processor.n_groups >= last_n_groups+10000:
        n_bam_aligns = merged_processor.stat_sum('n_bam_aligns')
        print >> sys.stderr, "Found %d reads across all BAM files (%d alignments processed -- %.2f per second)" % (merged_processor.n_groups, n_bam_aligns, (n_bam_aligns-last_n_bam_aligns)/bam_reading_stopwatch.reset())
        last_n_bam_aligns = n_bam_aligns
        last_n_groups = merged_processor.n_groups

BAMProcessor.n_paired_alignments_after_filter = property(lambda self: self.n_paired_alignments-self.n_non_exonic_bam_aligns-self.n_different_genes_pair-(self.n_partially_exonic if options.strict else 0))
attr_names = ['n_bam_aligns', 'n_singletons', 'n_multiple_hits', 'n_paired_alignments', 'n_non_exonic_bam_aligns', 'n_different_genes_pair', 'n_partially_exonic', 'n_fully_exonic', 'n_paired_alignments_after_filter', 'n_gene_match']
for stat_name in attr_names:
    setattr(merged_processor, stat_name, merged_processor.stat_sum(stat_name))

print >> sys.stderr, "Finished reading the BAM files in %.2f seconds" % main_stopwatch.reset()
print >> sys.stderr, "%d alignments in total across all %d BAM files" % (merged_processor.n_bam_aligns, len(bam_files))
print >> sys.stderr, "\t%d discarded (%.2f%%) - singletons " % (merged_processor.n_singletons, 100.*merged_processor.n_singletons/merged_processor.n_bam_aligns)
print >> sys.stderr, "\t%d discarded (%.2f%%) - multiple hits (NH != 1)" % (merged_processor.n_multiple_hits, 100.*merged_processor.n_multiple_hits/merged_processor.n_bam_aligns)
print >> sys.stderr, "%d paired alignments in total" % merged_processor.n_paired_alignments
if not options.no_gtf_filter:
    print >> sys.stderr, "\t%d discarded (%.2f%%) - both not exonic" % (merged_processor.n_non_exonic_bam_aligns, 100.*merged_processor.n_non_exonic_bam_aligns/merged_processor.n_paired_alignments)
    print >> sys.stderr, "\t%d discarded (%.2f%%) - multiple genes hit" % (merged_processor.n_different_genes_pair, 100.*merged_processor.n_different_genes_pair/merged_processor.n_paired_alignments)
    print >> sys.stderr, "\t%d %s (%.2f%%) - only one not exonic" % (merged_processor.n_partially_exonic, 'discarded' if options.strict else 'kept', 100.*merged_processor.n_partially_exonic/merged_processor.n_paired_alignments)
    print >> sys.stderr, "\t%d kept (%.2f%%) - both exonic in same gene" % (merged_processor.n_fully_exonic, 100.*merged_processor.n_fully_exonic/merged_processor.n_paired_alignments)
    print >> sys.stderr, "%d paired alignments in total after GTF filtering" % merged_processor.n_paired_alignments_after_filter
print >> sys.stderr, "%d unique read names" % merged_processor.n_groups
if merged_processor.n_groups and (not options.no_gtf_filter):
    print >> sys.stderr, "Gene name assignment statistics"
    print >> sys.stderr, "\t%d reads (%.2f%%): single candidate" % (merged_processor.n_unique_groups, 100.*merged_processor.n_unique_groups/merged_processor.n_groups)
    if merged_processor.n_best_groups:
        print >> sys.stderr, "\t%d reads (%.2f%%): multiple candidates, lowest NM score selected" % (merged_processor.n_best_groups, 100.*merged_processor.n_best_groups/merged_processor.n_groups)
    if merged_processor.n_ambiguous_groups:
        print >> sys.stderr, "\t%d reads (%.2f%%): multiple candidates, tie - %d candidates listed (%.2f per read on average)" % (merged_processor.n_ambiguous_groups, 100.*merged_processor.n_ambiguous_groups/merged_processor.n_groups, n_lines-merged_processor.n_unique_groups-merged_processor.n_best_groups, float(n_lines-merged_processor.n_unique_groups-merged_processor.n_best_groups)/merged_processor.n_ambiguous_groups)

if options.log_file:
    with open(options.log_file, 'w') as fh:
        print >> fh, '\t'.join(['# filename'] + attr_names)
        for (filename, bam_processor) in zip(bam_files, bam_processors):
            print >> fh, '\t'.join([filename] + [str(getattr(bam_processor, attr_name)) for attr_name in attr_names])


#Return a non-zero code if we couldn't find any groups
if not merged_processor.n_groups:
    sys.exit(1)
