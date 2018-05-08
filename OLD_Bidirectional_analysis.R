# Data reading
setwd("~/ProjectX")
matrix_plus <- read.table("EM.plus.tab", h = F, stringsAsFactors = F)
matrix_minus <- read.table("EM.minus.tab", h = F, stringsAsFactors = F)

# Column renaming
names(matrix_plus) <- c("POSITION", "SIZE", "C0", "C1", "C2", "C3", "KD0", "KD1", "KD2", "KD3")
names(matrix_minus) <- c("POSITION", "SIZE", "C0", "C1", "C2", "C3", "KD0", "KD1", "KD2", "KD3")

# Separate position into chromosome, start and end by pattern finding
myregex <- "([0-9A-z]*):([0-9]*)-([0-9]*)"

matrix_plus$CHR <- gsub(pattern = myregex, "\\1", x = matrix_plus$POSITION)
matrix_plus$START <- as.numeric(gsub(pattern = myregex, "\\2", x = matrix_plus$POSITION))
matrix_plus$END <- as.numeric(gsub(pattern = myregex, "\\3", x = matrix_plus$POSITION))

matrix_minus$CHR <- gsub(pattern = myregex, "\\1", x = matrix_minus$POSITION)
matrix_minus$START <- as.numeric(gsub(pattern = myregex, "\\2", x = matrix_minus$POSITION))
matrix_minus$END <- as.numeric(gsub(pattern = myregex, "\\3", x = matrix_minus$POSITION))

# Establish midpoints
matrix_plus$MIDS <- (matrix_plus$START + matrix_plus$END) / 2
matrix_minus$MIDS <- (matrix_minus$START + matrix_minus$END) / 2

# Correction factor based on window size (merging process originated windows longer than 400)
matrix_plus$CF <- 400 / matrix_plus$SIZE
matrix_minus$CF <- 400 / matrix_minus$SIZE

matrix_plus <- cbind(matrix_plus$POSITION, matrix_plus$MIDS, matrix_plus$CHR, 
                     matrix_plus[, c(3:10)] * matrix_plus$CF)
matrix_minus <- cbind(matrix_minus$POSITION, matrix_minus$MIDS, matrix_minus$CHR, 
                     matrix_minus[, c(3:10)] * matrix_minus$CF)

# Column renaming
names(matrix_plus) <- c("POSITION", "MIDPOINTS", "CHR", "C0", "C1", "C2", "C3", "KD0", "KD1", "KD2", "KD3")
names(matrix_minus) <- c("POSITION", "MIDPOINTS", "CHR", "C0", "C1", "C2", "C3", "KD0", "KD1", "KD2", "KD3")

# Chromosomes present in each expression matrix (the same in both):
# 2L, 2LHet, 2R, 2RHet, 3L, 3LHet, 3R, 3RHet, 4, U, Uextra, X, XHet
levels(matrix_plus$CHR); levels(matrix_minus$CHR)


# BOTH REPLICATES OF NF4T_GFP
NF4T_GFP_plus <- matrix_plus[, c(1:5)]
NF4T_GFP_minus <- matrix_minus[, c(1:5)]

# Both strands with a replicate expression level > 0 TPM
NF4T_GFP_plus <- NF4T_GFP_plus[apply(NF4T_GFP_plus[, c(4, 5)], MARGIN = 1, function(x) all(x > 0)), ]
NF4T_GFP_minus <- NF4T_GFP_minus[apply(NF4T_GFP_minus[, c(4, 5)], MARGIN = 1, function(x) all(x > 0)), ]

# Addition of average TPMs
NF4T_GFP_plus$PLUS <- rowMeans(NF4T_GFP_plus[, c(4, 5)])
NF4T_GFP_minus$MINUS <- rowMeans(NF4T_GFP_minus[, c(4, 5)])

# NF4T_GFP DHSs candidates of bidirectional transcription: 4948
NF4T_GFP_both <- merge(NF4T_GFP_plus[, c(1, 2, 3, 6)], NF4T_GFP_minus[, c(1, 6)], by = "POSITION") 

# Strand bias calculation
NF4T_GFP_both$BIAS <- NF4T_GFP_both$PLUS / (NF4T_GFP_both$PLUS + NF4T_GFP_both$MINUS)

# Classification of the directionality condition (U or B)
NF4T_GFP_both$TYPE <- ifelse(NF4T_GFP_both$BIAS <= 0.25 | NF4T_GFP_both$BIAS >= 0.75, "U", "B")

# Keep only bidirectionally transcribed DHS
NF4T_GFP_both <- NF4T_GFP_both[NF4T_GFP_both$TYPE == "B", ]


# BOTH REPLICATES OF NF4T_KD
NF4T_KD_plus <- matrix_plus[, c(1, 2, 3, 6, 7)]
NF4T_KD_minus <- matrix_minus[, c(1, 2, 3, 6, 7)]

# Both strands with a replicate expression level > 0 TPM
NF4T_KD_plus <- NF4T_KD_plus[apply(NF4T_KD_plus[, c(4, 5)], MARGIN = 1, function(x) all(x > 0)), ]
NF4T_KD_minus <- NF4T_KD_minus[apply(NF4T_KD_minus[, c(4, 5)], MARGIN = 1, function(x) all(x > 0)), ]

# Addition of average TPMs
NF4T_KD_plus$PLUS <- rowMeans(NF4T_KD_plus[, c(4, 5)])
NF4T_KD_minus$MINUS <- rowMeans(NF4T_KD_minus[, c(4, 5)])

# NF4T_KD DHSs candidates of bidirectional transcription: 7756
NF4T_KD_both <- merge(NF4T_KD_plus[, c(1, 2, 3, 6)], NF4T_KD_minus[, c(1, 6)], by = "POSITION") 

# Strand bias calculation
NF4T_KD_both$BIAS <- NF4T_KD_both$PLUS / (NF4T_KD_both$PLUS + NF4T_KD_both$MINUS)

# Classification of the directionality condition (U or B)
NF4T_KD_both$TYPE <- ifelse(NF4T_KD_both$BIAS <= 0.25 | NF4T_KD_both$BIAS >= 0.75, "U", "B")

# Keep only bidirectionally transcribed DHS
NF4T_KD_both <- NF4T_KD_both[NF4T_KD_both$TYPE == "B", ]


# BOTH REPLICATES OF QLAA_GFP
QLAA_GFP_plus <- matrix_plus[, c(1, 2, 3, 8, 9)]
QLAA_GFP_minus <- matrix_minus[, c(1, 2, 3, 8, 9)]

# Both strands with a replicate expression level > 0 TPM
QLAA_GFP_plus <- QLAA_GFP_plus[apply(QLAA_GFP_plus[, c(4, 5)], MARGIN = 1, function(x) all(x > 0)), ]
QLAA_GFP_minus <- QLAA_GFP_minus[apply(QLAA_GFP_minus[, c(4, 5)], MARGIN = 1, function(x) all(x > 0)), ]

# Addition of average TPMs
QLAA_GFP_plus$PLUS <- rowMeans(QLAA_GFP_plus[, c(4, 5)])
QLAA_GFP_minus$MINUS <- rowMeans(QLAA_GFP_minus[, c(4, 5)])

# QLAA_GFP DHSs candidates of bidirectional transcription: 5596
QLAA_GFP_both <- merge(QLAA_GFP_plus[, c(1, 2, 3, 6)], QLAA_GFP_minus[, c(1, 6)], by = "POSITION") 

# Strand bias calculation
QLAA_GFP_both$BIAS <- QLAA_GFP_both$PLUS / (QLAA_GFP_both$PLUS + QLAA_GFP_both$MINUS)

# Classification of the directionality condition (U or B)
QLAA_GFP_both$TYPE <- ifelse(QLAA_GFP_both$BIAS <= 0.25 | QLAA_GFP_both$BIAS >= 0.75, "U", "B")

# Keep only bidirectionally transcribed DHS
QLAA_GFP_both <- QLAA_GFP_both[QLAA_GFP_both$TYPE == "B", ]


# BOTH REPLICATES OF QLAA_KD
QLAA_KD_plus <- matrix_plus[, c(1, 2, 3, 10, 11)]
QLAA_KD_minus <- matrix_minus[, c(1, 2, 3, 10, 11)]

# Both strands with a replicate expression level > 0 TPM
QLAA_KD_plus <- QLAA_KD_plus[apply(QLAA_KD_plus[, c(4, 5)], MARGIN = 1, function(x) all(x > 0)), ]
QLAA_KD_minus <- QLAA_KD_minus[apply(QLAA_KD_minus[, c(4, 5)], MARGIN = 1, function(x) all(x > 0)), ]

# Addition of average TPMs
QLAA_KD_plus$PLUS <- rowMeans(QLAA_KD_plus[, c(4, 5)])
QLAA_KD_minus$MINUS <- rowMeans(QLAA_KD_minus[, c(4, 5)])

# QLAA_KD DHSs candidates of bidirectional transcription: 7494
QLAA_KD_both <- merge(QLAA_KD_plus[, c(1, 2, 3, 6)], QLAA_KD_minus[, c(1, 6)], by = "POSITION") 

# Strand bias calculation
QLAA_KD_both$BIAS <- QLAA_KD_both$PLUS / (QLAA_KD_both$PLUS + QLAA_KD_both$MINUS)

# Classification of the directionality condition (U or B)
QLAA_KD_both$TYPE <- ifelse(QLAA_KD_both$BIAS <= 0.25 | QLAA_KD_both$BIAS >= 0.75, "U", "B")

# Keep only bidirectionally transcribed DHS
QLAA_KD_both <- QLAA_KD_both[QLAA_KD_both$TYPE == "B", ]


## EXOSOME SENSITIVITY of divergent DHSs in NF4T condition

# Clear "midpoints" and "chr" from NF4T-KD-BOTH DF
NF4T_KD_both_red <- NF4T_KD_both[, c(1, 4, 5, 6, 7)]

# Merge reduced DF including P and M parameters before and after exosome KD: 3353 (some positions lost)
NF4T_exosome <- merge(NF4T_GFP_both, NF4T_KD_both_red, by = "POSITION", suffixes = c(".GFP", ".KD"))

# Add the P_sensitivity and M_sensitivity columns
NF4T_exosome$P_sensitivity <- (NF4T_exosome$PLUS.KD - NF4T_exosome$PLUS.GFP) / NF4T_exosome$PLUS.KD
NF4T_exosome$M_sensitivity <- (NF4T_exosome$MINUS.KD - NF4T_exosome$MINUS.GFP) / NF4T_exosome$MINUS.KD

# Classification of the exosome sensitivity condition: stable (Stb), unstable (Uns) or none (Unclassified)
NF4T_exosome$P_type <- ifelse(NF4T_exosome$P_sensitivity <= 0.25, "Stb",
                              ifelse(NF4T_exosome$P_sensitivity >= 0.75, "Uns", "Unclassified"))
NF4T_exosome$M_type <- ifelse(NF4T_exosome$M_sensitivity <= 0.25, "Stb",
                              ifelse(NF4T_exosome$M_sensitivity >= 0.75, "Uns", "Unclassified"))


## EXOSOME SENSITIVITY of divergent DHSs in QLAA condition

# Clear "midpoints" and "chr" from QLAA KD DF
QLAA_KD_both_red <- QLAA_KD_both[, c(1, 4, 5, 6, 7)]

# Merge reduced DF including P and M parameters before and after exosome KD: 5095 (some positions lost)
QLAA_exosome <- merge(QLAA_GFP_both, QLAA_KD_both_red, by = "POSITION", suffixes = c(".GFP", ".KD"))

# Add the P_sensitivity and M_sensitivity columns
QLAA_exosome$P_sensitivity <- (QLAA_exosome$PLUS.KD - QLAA_exosome$PLUS.GFP) / QLAA_exosome$PLUS.KD
QLAA_exosome$M_sensitivity <- (QLAA_exosome$MINUS.KD - QLAA_exosome$MINUS.GFP) / QLAA_exosome$MINUS.KD

# Classification of the exosome sensitivity condition: stable (Stb), unstable (Uns) or none (Unclassified)
QLAA_exosome$P_type <- ifelse(QLAA_exosome$P_sensitivity <= 0.25, "Stb",
                              ifelse(QLAA_exosome$P_sensitivity >= 0.75, "Uns", "Unclassified"))
QLAA_exosome$M_type <- ifelse(QLAA_exosome$M_sensitivity <= 0.25, "Stb",
                              ifelse(QLAA_exosome$M_sensitivity >= 0.75, "Uns", "Unclassified"))


## Export complete DF as text files
write.table(NF4T_exosome, file = "NF4T_BD.txt", row.names = F, quote = F)
write.table(QLAA_exosome, file = "QLAA_BD.txt", row.names = F, quote = F)


                              