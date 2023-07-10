
hap_drb_name <- c(rep("DR1", 2), 
             rep("DR51", 2), 
             rep("DR52", 5), 
             rep("DR53", 3), 
             rep("DR8", 1) )
DRB1 <- c("01", "10", 
          "15", "16",
          "03", "11", "12", "13", "14",
          "04", '07', "09",
          "08")
hap_drb <- c(rep("DRB6;DRB9", 2), 
             rep("DRB5;DRB6;DRB9", 2), 
             rep("DRB2;DRB3;DRB9", 5), 
             rep("DRB4;DRB7;DRB8;DRB9", 3), 
             rep("DRB3;DRB9", 1) )

drb_haplotypes <- data.frame(Hap_name=hap_drb_name,
                             DRB1=DRB1,
                             DRB_par=hap_drb)
