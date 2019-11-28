library(tidyverse)
library(treeWAS)
library(ape)


meta_raw <- read_csv("mtb_meta.csv")
isos_raw <- read_table("mtb_iso_list.txt")

meta <- select(meta_raw, isolates, INH, RIF, EMB, STR, CIP, CYS, CAP, ETA, KAN, OFLX, PAS, PZA, AMK, MOXI, PRO, CLO, RBU, LEVO)
meta[meta == "NULL"] <- NA

isos_meta <- semi_join(meta, isos_raw, by = "isolates")

mtb_phen_inh <- as_vector(isos_meta$INH)
names(mtb_phen_inh) <- isos_meta$isolates

mtb_phen_rif <- as_vector(isos_meta$RIF)
names(mtb_phen_rif) <- isos_meta$isolates

mtb_phen_emb <- as_vector(isos_meta$EMB)
names(mtb_phen_emb) <- isos_meta$isolates

mtb_phen_str <- as_vector(isos_meta$STR)
names(mtb_phen_str) <- isos_meta$isolates

mtb_phen_cip <- as_vector(isos_meta$CIP)
names(mtb_phen_cip) <- isos_meta$isolates

mtb_phen_cys <- as_vector(isos_meta$CYS)
names(mtb_phen_cys) <- isos_meta$isolates

mtb_phen_cap <- as_vector(isos_meta$CAP)
names(mtb_phen_cap) <- isos_meta$isolates

mtb_phen_eta <- as_vector(isos_meta$ETA)
names(mtb_phen_eta) <- isos_meta$isolates

mtb_phen_kan <- as_vector(isos_meta$KAN)
names(mtb_phen_kan) <- isos_meta$isolates

mtb_phen_oflx <- as_vector(isos_meta$OFLX)
names(mtb_phen_oflx) <- isos_meta$isolates

mtb_phen_pas <- as_vector(isos_meta$PAS)
names(mtb_phen_pas) <- isos_meta$isolates

mtb_phen_pza <- as_vector(isos_meta$PZA)
names(mtb_phen_pza) <- isos_meta$isolates

mtb_phen_amk <- as_vector(isos_meta$AMK)
names(mtb_phen_amk) <- isos_meta$isolates

mtb_phen_moxi <- as_vector(isos_meta$MOXI)
names(mtb_phen_moxi) <- isos_meta$isolates

mtb_phen_pro <- as_vector(isos_meta$PRO)
names(mtb_phen_pro) <- isos_meta$isolates

mtb_phen_clo <- as_vector(isos_meta$CLO)
names(mtb_phen_clo) <- isos_meta$isolates

mtb_phen_rbu <- as_vector(isos_meta$RBU)
names(mtb_phen_rbu) <- isos_meta$isolates

mtb_phen_levo <- as_vector(isos_meta$LEVO)
names(mtb_phen_levo) <- isos_meta$isolates

mtb_dna <- read.dna(file = "mtb_snv_align.fasta", format = "fasta")
mtb_snps <- DNAbin2genind(mtb_dna)@tab
mtb_tree <- read.tree(file = "mtb_tree.nhx")

inh <- treeWAS(snps = mtb_snps,
               phen = mtb_phen_inh,
               tree = mtb_tree,
               seed = 1)

rif <- treeWAS(snps = mtb_snps,
                   phen = rif,
                   tree = mtb_tree,
                   seed = 1)

emb <- treeWAS(snps = mtb_snps,
               phen = mtb_phen_emb,
               tree = mtb_tree,
               seed = 1)

str <- treeWAS(snps = mtb_snps,
               phen = mtb_phen_str,
               tree = mtb_tree,
               seed = 1)

cip  <- treeWAS(snps = mtb_snps,
                   phen = mtb_phen_cip,
                   tree = mtb_tree,
                   seed = 1)

cys  <- treeWAS(snps = mtb_snps,
                   phen = mtb_phen_cys,
                   tree = mtb_tree,
                   seed = 1)

cap  <- treeWAS(snps = mtb_snps,
                phen = mtb_phen_cap,
                tree = mtb_tree,
                seed = 1)

eta  <- treeWAS(snps = mtb_snps,
                phen = mtb_phen_eta,
                tree = mtb_tree,
                seed = 1)

kan  <- treeWAS(snps = mtb_snps,
                phen = mtb_phen_kan,
                tree = mtb_tree,
                seed = 1)

oflx   <- treeWAS(snps = mtb_snps,
                  phen = mtb_phen_oflx,
                  tree = mtb_tree,
                  seed = 1)

pas  <- treeWAS(snps = mtb_snps,
                phen = mtb_phen_pas,
                tree = mtb_tree,
                seed = 1)

pza  <- treeWAS(snps = mtb_snps,
                phen = mtb_phen_pza,
                tree = mtb_tree,
                seed = 1)

amk  <- treeWAS(snps = mtb_snps,
                phen = mtb_phen_amk,
                tree = mtb_tree,
                seed = 1)

moxi  <- treeWAS(snps = mtb_snps,
                 phen = mtb_phen_moxi,
                 tree = mtb_tree,
                 seed = 1)

pro  <- treeWAS(snps = mtb_snps,
                phen = mtb_phen_pro,
                tree = mtb_tree,
                seed = 1)

clo  <- treeWAS(snps = mtb_snps,
                phen = mtb_phen_clo,
                tree = mtb_tree,
                seed = 1)

rbu  <- treeWAS(snps = mtb_snps,
                phen = mtb_phen_rbu,
                tree = mtb_tree,
                seed = 1)

levo  <- treeWAS(snps = mtb_snps,
                 phen = mtb_phen_levo,
                 tree = mtb_tree,
                 seed = 1)
