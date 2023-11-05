## code to prepare `DATASET` dataset goes here

usethis::use_data(DATASET, overwrite = TRUE)


# Internal data to be used

# Preparing human transcription factor data

humanTF <- read.csv("./data-raw/DatabaseExtract_v_1.01.csv",
  header = TRUE,
  row.names = 1
)
humanTF <- humanTF[humanTF$Is.TF. == "Yes", ]
humanTF <- humanTF %>% dplyr::select(
  HGNC.symbol, Ensembl.ID, EntrezGene.ID,
  EntrezGene.Description, DBD, PDB.ID
)

rownames(humanTF) <- humanTF$HGNC.symbol

TFHsapiens <- humanTF


# Document human small RNAs
write.csv(miRNAHsapiens, "./data-raw/miRNAHsapiens_COMPSRA.csv", row.names = FALSE)
write.csv(snoRNAHsapiens, "./data-raw/snoRNAHsapiens_COMPSRA.csv", row.names = FALSE)
write.csv(snRNAHsapiens, "./data-raw/snRNAHsapiens_COMPSRA.csv", row.names = FALSE)
write.csv(tRNAHsapiens, "./data-raw/tRNAHsapiens_COMPSRA.csv", row.names = FALSE)
write.csv(piRNAHsapiens, "./data-raw/piRNAHsapiens_COMPSRA.csv", row.names = FALSE)
write.csv(circRNAHsapiens, "./data-raw/circRNAHsapiens_COMPSRA.csv", row.names = FALSE)




# Saving internal data
usethis::use_data(TFHsapiens, miRNAHsapiens, snoRNAHsapiens, snRNAHsapiens,
  tRNAHsapiens, piRNAHsapiens, circRNAHsapiens,
  internal = TRUE, overwrite = TRUE
)
