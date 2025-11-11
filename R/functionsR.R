
# Functions ---------------------------------------------------------------

# Recognize user defined functions (f.abc) and arguments (a.xyz)
# Plot functions are at the script end (f.p.scatter) together with
# the respective plot export function.



# Import DIA-NN matrix from specified folder into data frame.
# Uses partial filename match, e.g."pg_matrix", "pr_matrix".
# requires excess matrices to be moved away. (Applies when using MBR)
# e.g. Prot0 <- f.import.matrix(folder_input, "pg_matrix", "\t")
f.import.matrix <- function(a.folder,
                            a.matrix,
                            a.sep) {
    read.csv(sep = a.sep,
             check.names = FALSE,
             list.files(a.folder,
                        pattern = a.matrix,
                        full.names = TRUE))
}




# In DIA-NN matrices, rename columns containing quantitative values
# from filepath to short (stripped) replicate name.
# Remove filepath and suffix, remainder is used as replicate name.
# e.g. Prot1 <- f.strip.sample.column.name(Prot1, col_ctr_Prot)
f.strip.sample.column.name <- function(a.df, a.columns) {

    names(a.df)[a.columns] <- gsub(".raw$", "", names(a.df)[a.columns])
    names(a.df)[a.columns] <- gsub(".dia$", "", names(a.df)[a.columns])
    names(a.df)[a.columns] <- gsub(".mzML$", "", names(a.df)[a.columns])
    names(a.df)[a.columns] <- gsub(".wiff$", "", names(a.df)[a.columns])
    names(a.df)[a.columns] <- gsub(".wiff2$", "", names(a.df)[a.columns])
    names(a.df)[a.columns] <- gsub(".d$", "", names(a.df)[a.columns])

    names(a.df)[a.columns] <- gsub(".raw$", "", names(a.df)[a.columns])
    names(a.df)[a.columns] <- gsub(".dia$", "", names(a.df)[a.columns])
    names(a.df)[a.columns] <- gsub(".mzML$", "", names(a.df)[a.columns])
    names(a.df)[a.columns] <- gsub(".wiff$", "", names(a.df)[a.columns])
    names(a.df)[a.columns] <- gsub(".wiff2$", "", names(a.df)[a.columns])
    names(a.df)[a.columns] <- gsub(".d$", "", names(a.df)[a.columns])

    names(a.df)[a.columns] <- gsub(".raw$", "", names(a.df)[a.columns])
    names(a.df)[a.columns] <- gsub(".dia$", "", names(a.df)[a.columns])
    names(a.df)[a.columns] <- gsub(".mzML$", "", names(a.df)[a.columns])
    names(a.df)[a.columns] <- gsub(".wiff$", "", names(a.df)[a.columns])
    names(a.df)[a.columns] <- gsub(".wiff2$", "", names(a.df)[a.columns])
    names(a.df)[a.columns] <- gsub(".d$", "", names(a.df)[a.columns])

    names(a.df)[a.columns] <- basename(names(a.df)[a.columns])
    # names(a.df)[a.columns] <- sub(".[^.]+$", "", colnames(a.df)[a.columns])
    return(a.df)
}




# Remove entries of MaxQuant contaminant fasta from an e.g. DIA-NN pg_matrix.
# e.g. Prot1 <- f.remove.MQ.cont(Prot1)
f.remove.MQ.cont <- function(a.df) {
    a.df <- a.df[!grepl("TREMBL", a.df[, "Genes"]),]
    a.df <- a.df[!grepl("SWISS-PROT", a.df[, "Genes"]),]
    a.df <- a.df[!grepl("REFSEQ", a.df[, "Protein.Ids"]),]
    a.df <- a.df[!grepl("ENSEMBL", a.df[, "Protein.Ids"]),]
    a.df <- a.df[!grepl("H-INV:HIT", a.df[, "Protein.Ids"]),]
    a.df <- a.df[!grepl("Streptavidin", a.df[, "Protein.Ids"]),]
}




# For each entry and condition, calculate stats of the replicate measurements
# such as count of non-0 and non-NA values, mean, standard deviation,
# coefficient of variation (CV - standard deviation in % of the mean)
# e.g. Prot1 <- f.row.stats(Prot1, col_ctr_Prot, col_exp_Prot)
f.row.stats <- function(a.df, a.col_ctr, a.col_exp) {
    a.df[, "ctr_count"] <- apply(!is.na(a.df[, a.col_ctr]), 1, sum)
    a.df[, "ctr_mean"]  <- rowMeans(a.df[, a.col_ctr], na.rm = TRUE)
    a.df[, "ctr_SD"]    <-
        rowSds(as.matrix(a.df[, a.col_ctr]), na.rm = TRUE)
    a.df[, "ctr_CV"]    <- a.df[, "ctr_SD"] / a.df[, "ctr_mean"] * 100

    a.df[, "exp_count"] <- apply(!is.na (a.df[, a.col_exp]), 1, sum)
    a.df[, "exp_mean"]  <- rowMeans(a.df[, a.col_exp], na.rm = TRUE)
    a.df[, "exp_SD"]    <-
        rowSds(as.matrix(a.df[, a.col_exp]), na.rm = TRUE)
    a.df[, "exp_CV"]    <- a.df[, "exp_SD"] / a.df[, "exp_mean"] * 100

    return(a.df)
}




# For LFQbenchmark experiments, perform a row of basic data sanitation
# to prepare for later following data filtering.
# This contains annotating entries with respective species
# and retaining only single-species entries.
# e.g. Prot1 <- f.LFQbenchmark.basics(Prot0, col_ctr_Prot, col_exp_Prot)
f.LFQbenchmark.basics <- function(a.df,
                                  a.col_ctr,
                                  a.col_exp) {


    # Replace 0 with NA. "0" is sometimes used as mock value.
    # The mock value should not enter any calculation
    # and the event should be seen as missing value.(NA)
    a.df[a.col_ctr][a.df[a.col_ctr] == 0] <- NA
    a.df[a.col_exp][a.df[a.col_exp] == 0] <- NA



    # Calculate log2FC,
    # can be slightly different (5^-12 units) from limma,
    # but is easier to reproduce.
    a.df[, "log2FC"] <- log2(a.df[, "exp_mean"] / a.df[, "ctr_mean"])



    # In LFQbenchmark experiments, assign the respective species
    # to species-specific entries and remove entries matching
    # multiple species or no species. Should also work on 2-species.
    a.df[, "Species"] <- NA


    a.df[, "Species"][grepl('_HUMAN', a.df[, "Protein.Names"]) &
                          !grepl('_ECOLI', a.df[, "Protein.Names"]) &
                          !grepl('_YEAST', a.df[, "Protein.Names"]) &
                          !grepl('_CAEEL', a.df[, "Protein.Names"])] <- 'Human'

    a.df[, "Species"][grepl('_YEAST', a.df[, "Protein.Names"]) &
                          !grepl('_HUMAN', a.df[, "Protein.Names"]) &
                          !grepl('_ECOLI', a.df[, "Protein.Names"]) &
                          !grepl('_CAEEL', a.df[, "Protein.Names"])] <- 'Yeast'

    a.df[, "Species"][grepl('_ECOLI', a.df[, "Protein.Names"]) &
                          !grepl('_HUMAN', a.df[, "Protein.Names"]) &
                          !grepl('_YEAST', a.df[, "Protein.Names"]) &
                          !grepl('_CAEEL', a.df[, "Protein.Names"])] <- 'E.coli'

    a.df[, "Species"][grepl('_CAEEL', a.df[, "Protein.Names"]) &
                          !grepl('_HUMAN', a.df[, "Protein.Names"]) &
                          !grepl('_YEAST', a.df[, "Protein.Names"]) &
                          !grepl('_ECOLI', a.df[, "Protein.Names"])] <- 'C.elegans'


    a.df <- a.df[complete.cases(a.df[, "Species"]), ]





    # Add column listing expected log2FC row-wise
    # to simplify calculations of summary_stats around the log2FC.

    a.df[, "expected_log2FC"] <- NA


    a.df[a.df[, "Species"]== "Human", "expected_log2FC"]      <- expFC_human
    a.df[a.df[, "Species"]== "Yeast", "expected_log2FC"]      <- expFC_yeast
    a.df[a.df[, "Species"]== "E.coli", "expected_log2FC"]     <- expFC_ecoli
    a.df[a.df[, "Species"]== "C.elegans", "expected_log2FC"]  <- expFC_celegans


    return(a.df)

}




# Acquire adjusted p-values for differential expression analysis
# for the protein group matrix with limma. Not executed on Precursor level.
# e.g. Prot3 <- f.limma.p_adj(Prot3, col_ctr_Prot, col_exp_Prot)
f.limma.p_adj <- function(a.df,
                          a.col_ctr,
                          a.col_exp) {

    # Limma to get BH adjusted p-values (p_adj) for differential expression.
    dcontra = matrix(c('exp', 'con'),
                     ncol = 2,
                     dimnames = list(NULL, c('CONDITION', 'CONDITION_REF'))) # defines in each rwo conditions that should be compared
    assay = as.matrix(a.df[, c(a.col_exp, a.col_ctr)]) # features x samples: Protein.Group x exp and con columns
    rownames(assay) = a.df[, 'Protein.Group']


    ## IF e.g. grps=c('exp','exp','exp','ctr','ctr','ctr')
    ## then for paired analysis use reps=c(1,2,3,1,2,3)
    ## Otherwise "reps" and "limma.rep" lines should be excluded from running.

    grps = c(rep('exp', length(a.col_exp)), rep('con', length(a.col_ctr))) # conditions of samples in assay


    # in case of paired analysis e.g. reps=c(1,2,3,1,2,3)
    # reps = c(1:length(a.col_exp), 1:length(a.col_ctr))


    # Limma expects logarithmized input intensities
    assay = log2(assay)

    limma_data = ExpressionSet(assay
                               ,
                               phenoData = AnnotatedDataFrame(data.frame(
                                   CONDITION = grps
                                   , row.names = colnames(assay)
                               ))
                               ,
                               featureData = AnnotatedDataFrame(
                                   data.frame(
                                       PROTEIN_IDS = a.df[, 'Protein.Ids']
                                       ,
                                       PROTEIN_NAMES = a.df[, 'Protein.Names']
                                       ,
                                       GENES = a.df[, 'Genes']
                                       ,
                                       row.names = rownames(assay)
                                   )
                               ))
    comparisons <-
        c('exp - con') # Conditions that should be compared, in the trivial case it's exp vs con
    limma.cond = factor(grps)

    ## in case of paired analysis enable limma.rep
    #limma.rep <- factor(reps)
    contrast.matrix <- model.matrix(~  0 + limma.cond
                                    # + limma.rep
    )
    colnames(contrast.matrix) <-
        gsub("limma.cond", "", colnames(contrast.matrix))

    limma.object <- eBayes(contrasts.fit(
        lmFit(limma_data, design = contrast.matrix)
        ,
        makeContrasts(contrasts = comparisons, levels = contrast.matrix)
    )
    ,
    trend = TRUE,
    robust = TRUE)
    res_limma = limma::topTable(
        limma.object,
        coef = comparisons[1],
        number = Inf,
        sort.by = 'none',
        adjust.method = "BH",
        confint = TRUE
    )


    # Add p_adj (limma result) to the main protein group data frame.
    a.df[, "p_adj"] <- res_limma[, 'adj.P.Val']
    a.df[, "p_adj"] <- as.numeric(format(a.df[, "p_adj"],
                                         scientific = FALSE,
                                         justified = "none"))
    return(a.df)


    # # For use outside of a function,
    # # clean up environment from unused limma elements.
    # rm(assay)
    # rm(contrast.matrix)
    # rm(dcontra)
    # rm(limma_data)
    # rm(limma.object)
    # rm(comparisons)
    # rm(grps)
    # rm(limma.cond)
    # rm(limma.rep)
    # rm(reps)
    # rm(res_limma)
}




# "Interpret" the p_adj from limma differential expression analysis,
# this generates the "confusion matrix" summary stats.
# !!! Hard-coded to the orientation of the example in the introduction.
# Please see introduction.
# e.g. Prot3 <- f.LFQbenchmark.limma.interpretation(Prot3, alpha_limma, limit_FC)
f.LFQbenchmark.limma.interpretation <-
    function(a.df, a.alpha_limma, a.limit_FC) {

        a.df[(a.df[, "Species"] == "Human"
              & (a.df[, "p_adj"] >= a.alpha_limma
                 | abs(a.df[, "log2FC"]) <= a.limit_FC)), "DE_result"] <-
            "true negative"

        a.df[(a.df[, "Species"] == "Human"
              & a.df[, "p_adj"] < a.alpha_limma
              & abs(a.df[, "log2FC"]) > a.limit_FC), "DE_result"] <-
            "false positive"



        a.df[(a.df[, "Species"] == "Yeast"
              & a.df[, "log2FC"] > +a.limit_FC
              & a.df[, "p_adj"] < a.alpha_limma), "DE_result"] <-
            "true positive"

        a.df[(a.df[, "Species"] == "Yeast"
              & a.df[, "log2FC"] < -a.limit_FC
              & a.df[, "p_adj"] < a.alpha_limma), "DE_result"] <-
            "false positive"

        a.df[(a.df[, "Species"] == "Yeast"
              & (a.df[, "p_adj"] >= a.alpha_limma
                 | abs(a.df[, "log2FC"]) <= a.limit_FC)), "DE_result"] <-
            "false negative"



        a.df[(a.df[, "Species"] == "E.coli"
              & a.df[, "log2FC"] < -a.limit_FC
              & a.df[, "p_adj"] < a.alpha_limma), "DE_result"] <-
            "true positive"

        a.df[(a.df[, "Species"] == "E.coli"
              & a.df[, "log2FC"] > +a.limit_FC
              & a.df[, "p_adj"] < a.alpha_limma), "DE_result"] <-
            "false positive"

        a.df[(a.df[, "Species"] == "E.coli"
              & (a.df[, "p_adj"] >= a.alpha_limma
                 | abs(a.df[, "log2FC"]) <= a.limit_FC)), "DE_result"] <-
            "false negative"



        a.df[(a.df[, "Species"] == "C.elegans"
              & a.df[, "log2FC"] < -a.limit_FC
              & a.df[, "p_adj"] < a.alpha_limma), "DE_result"] <-
            "true positive"

        a.df[(a.df[, "Species"] == "C.elegans"
              & a.df[, "log2FC"] > +a.limit_FC
              & a.df[, "p_adj"] < a.alpha_limma), "DE_result"] <-
            "false positive"

        a.df[(a.df[, "Species"] == "C.elegans"
              & (a.df[, "p_adj"] >= a.alpha_limma
                 | abs(a.df[, "log2FC"]) <= a.limit_FC)), "DE_result"] <-
            "false negative"


        return(a.df)

    }




# Tailing and asymmetry factor for up- and down-regulated proteins (according to species mixtures)
# Different calculations to acquire tailing towards the outside away from log2FC of 0 (or vice versa),
# in both directions.
# E.g. f.asym.down-regulated(Prot3, "E.coli")
f.asym.downregulated <- function(a.df,
                                 a.species) {

    # if statement required to prevent errors when handling 2-species benchmarks.
    if(length(a.df[a.df[, "Species"] == a.species, "log2FC"]) > 2) {

        # 1-dimensional data like log2 fold-changes
        # to data frames containing x and y coordinates of density plot function
        data <- a.df[a.df[, "Species"] == a.species, "log2FC"]
        density <- density(data)
        density <- as.data.frame(cbind(density$x, density$y))
        colnames(density) <- c("x", "y")

        # ymax and respective x coordinate called C
        ymax <- density$y[which.max(density$y)]
        C <- density$x[which.max(density$y)]

        # subset density coordinates by lefthand and righthand of modal
        # to differentiate multiple x coordinates matching one y coordinate.
        density_left <- density[density$x < C,]
        density_right <- density[density$x > C,]

        # A and B are the index, x, and y coordinates at a given height relative to ymax.
        A1 <- density_left[which.min(abs(density_left$y - ymax*0.10)),]
        B1 <- density_right[which.min(abs(density_right$y - ymax*0.10)),]

        A2 <- density_left[which.min(abs(density_left$y - ymax*0.05)),]
        B2 <- density_right[which.min(abs(density_right$y - ymax*0.05)),]

        # Factors calculated directional towards the outside away from log2FC of 0.
        Asymmetry_Factor <- (A1$x - C) / (C - B1$x)
        Tailing_Factor <- (A2$x - B2$x) / (2* (C - B2$x))

        return(cbind(Tailing_Factor, Asymmetry_Factor))

    } else {

        return(cbind(NA, NA))}

}

# Asymmetry and tailing calculated for up-regulated log2 fold-changes.
# # E.g. f.asym.up-regulated(Prot3, "Yeast")
f.asym.upregulated <- function(a.df,
                               a.species) {

    # if statement required to prevent errors when handling 2-species benchmarks.
    if(length(a.df[a.df[, "Species"] == a.species, "log2FC"]) > 2) {

        # 1-dimensional data like log2 fold-changes
        # to data frames containing x and y coordinates of density plot function
        data <- a.df[a.df[, "Species"] == a.species, "log2FC"]
        density <- density(data)
        density <- as.data.frame(cbind(density$x, density$y))
        colnames(density) <- c("x", "y")

        # ymax and respective x coordinate called C
        ymax <- density$y[which.max(density$y)]
        C <- density$x[which.max(density$y)]

        # subset density coordinates by leftside and rightside of modal
        # to differentiate x coordinates matchign one y coordiante
        density_left <- density[density$x < C,]
        density_right <- density[density$x > C,]

        # A and B are the index, x, and y coordinates at a given height relative to ymax.
        A1 <- density_left[which.min(abs(density_left$y - ymax*0.10)),]
        B1 <- density_right[which.min(abs(density_right$y - ymax*0.10)),]

        A2 <- density_left[which.min(abs(density_left$y - ymax*0.05)),]
        B2 <- density_right[which.min(abs(density_right$y - ymax*0.05)),]

        # Factors calculated directional towards the outside away from log2FC of 0.
        Asymmetry_Factor <- (B1$x - C) / (C - A1$x)
        Tailing_Factor <- (B2$x - A2$x) / (2* (C - A2$x))

        return(cbind(Tailing_Factor, Asymmetry_Factor))

    } else {

        return(cbind(NA, NA))}

}




# For processed LFQbenchmark data, collect main summary statistics
# into a new data frame for export and plot subtitles.
# They describe sensitivity, precision, accuracy, reliability, etc..
# Next to plots they are the core result of this script.
# Mainly use true positives (TP), deFDR, and Asymmetry_Factor to evaluate the performance.
# See introduction for interpretation.
f.LFQbenchmark.summary.stats <-
    function(a.df01, a.df02, a.df03, a.df04) {
        summary_stats <- data.frame(matrix(ncol = 0, nrow = 1))


        # R script filter variables to avoid confusion
        # when collecting stats of multiple benchmarks in excel.
        summary_stats[, "Variables"] <-
            paste0(script_version,",",
                   rep_min, "of", rep_max,",",
                   limit_CV,",",
                   limit_FC,",",
                   alpha_limma)


        # deFDR (differential expression FDR) (part of confusion matrix)
        summary_stats[, "deFDR"] <-
            round(digits = 2,
                  100 * length(which(a.df02[, "DE_result"] == "false positive")) /
                      (length(which(
                          a.df02[, "DE_result"] == "false positive"
                      ))
                      + length(which(
                          a.df02[, "DE_result"] == "true positive"
                      ))))


        # true positive count
        summary_stats[, "TP"] <-
            length(which(a.df02[, "DE_result"] == "true positive"))


        # Number of confidently detected and quantified protein groups,
        # after filtering for missingness and missingness + CV.
        summary_stats[, "Prot_ID"] <- nrow(a.df01)
        summary_stats[, "Prot_Quant"] <- nrow(a.df02)


        # Asymmetry_Factor to evaluate ratio compression or extension.
        summary_stats[, "Prot_Asymmetry_E.coli"] <-
            round( digits =2,
                   as.numeric(subset(asymmetry,
                                     Species == "E.coli" &
                                         Group == "protein group",
                                     select = "Asymmetry_Factor")))

        summary_stats[, "Prot_Asymmetry_Yeast"] <-
            round( digits =2,
                   as.numeric(subset(asymmetry,
                                     Species == "Yeast" &
                                         Group == "protein group",
                                     select = "Asymmetry_Factor")))


        # Coefficient of variation to describe precision.
        # For ca. 3 replicates aim at average and median at or below 5%
        # after filtering for <20%.
        summary_stats[, "Prot_CV_Mean"]   <-
            round(digits = 2,  mean(c(a.df02[, "ctr_CV"], a.df02[, "exp_CV"])))
        summary_stats[, "Prot_CV_Median"] <-
            round(digits = 2, median(c(a.df02[, "ctr_CV"], a.df02[, "exp_CV"])))


        # Accuracy is average distance of measured to expected log2FC.
        summary_stats[, "Prot_Accuracy"] <-
            round(digits = 2,
                  mean(abs(a.df02[, "log2FC"] - a.df02[, "expected_log2FC"])))


        # Trueness is distance between measurement medians
        # and respective expected fold-changes.
        # (cumulative as shifts are important to spot)
        # Several conditions can lead to a high value
        # e.g. ratio compression or erroneous normalization or even both together.
        summary_stats[, "Prot_Trueness"] <-
            round(
                digits = 2,
                sum(na.rm = TRUE,
                    c(abs(median(a.df02[a.df02[, "Species"] == "Human"  , "log2FC"]) - expFC_human),
                      abs(median(a.df02[a.df02[, "Species"] == "Yeast"  , "log2FC"]) - expFC_yeast),
                      abs(median(a.df02[a.df02[, "Species"] == "E.coli" , "log2FC"]) - expFC_ecoli),
                      abs(median(a.df02[a.df02[, "Species"] == "C.elegans"  , "log2FC"]) - expFC_celegans)
                    )))


        # Dispersion is the average distance of individual measurements around the respective median.
        summary_stats[, "Prot_Dispersion"] <-
            round(digits = 2,
                  mean(na.rm = TRUE,
                       c(
                           abs(a.df02[a.df02[, "Species"] == "Human"   , "log2FC"] - median(a.df02[a.df02[, "Species"] == "Human"  , "log2FC"])),
                           abs(a.df02[a.df02[, "Species"] == "Yeast"   , "log2FC"] - median(a.df02[a.df02[, "Species"] == "Yeast"  , "log2FC"])),
                           abs(a.df02[a.df02[, "Species"] == "E.coli"  , "log2FC"] - median(a.df02[a.df02[, "Species"] == "E.coli" , "log2FC"])),
                           abs(a.df02[a.df02[, "Species"] == "C.elegans"  , "log2FC"] - median(a.df02[a.df02[, "Species"] == "C.elegans" , "log2FC"]))
                       )))


        # Other than TP, continued confusion matrix of
        # differential expression result interpretation.
        summary_stats[, "FP"] <-
            length(which(a.df02[, "DE_result"] == "false positive"))
        summary_stats[, "TN"] <-
            length(which(a.df02[, "DE_result"] == "true negative"))
        summary_stats[, "FN"] <-
            length(which(a.df02[, "DE_result"] == "false negative"))


        summary_stats[, "Sensitivity"] <-
            round(digits = 2, 100 * length(which(a.df02[, "DE_result"] == "true positive")) /
                      (length(which(
                          a.df02[, "DE_result"] == "true positive"
                      )) + length(which(
                          a.df02[, "DE_result"] == "false negative"
                      ))))

        summary_stats[, "Specificity"] <-
            round(digits = 2, 100 * length(which(a.df02[, "DE_result"] == "true negative")) /
                      (length(which(
                          a.df02[, "DE_result"] == "true negative"
                      )) + length(which(
                          a.df02[, "DE_result"] == "false positive"
                      ))))






        # Stats partially repeated on precursor level as for protein groups above.
        summary_stats[, "Prec_ID"] <- nrow(a.df03)
        summary_stats[, "Prec_Quant"] <- nrow(a.df04)


        # Asymmetry_Factor to evaluate ratio compression or extension.
        summary_stats[, "Prec_Asymmetry_E.coli"] <-
            round( digits =2,
                   as.numeric(subset(asymmetry,
                                     Species == "E.coli" &
                                         Group == "precursor",
                                     select = "Asymmetry_Factor")))

        summary_stats[, "Prec_Asymmetry_Yeast"] <-
            round( digits =2,
                   as.numeric(subset(asymmetry,
                                     Species == "Yeast" &
                                         Group == "precursor",
                                     select = "Asymmetry_Factor")))



        # Precursor-level precision by coefficient of variation.
        summary_stats[, "Prec_CV_Mean"]   <-  round(digits = 2,
                                                    mean(c(a.df04[, "ctr_CV"], a.df04[, "exp_CV"])))
        summary_stats[, "Prec_CV_Median"] <-  round(digits = 2,
                                                    median(c(a.df04[, "ctr_CV"], a.df04[, "exp_CV"])))


        # Accuracy is average distance of measured to expected log2FC.
        summary_stats[, "Prec_Accuracy"] <-
            round(digits = 2,
                  mean(abs(a.df04[, "log2FC"] - a.df04[, "expected_log2FC"])))


        # Trueness is distance between measurement medians
        # and respective expected fold-changes (cumulative as shifts are important to spot)
        summary_stats[, "Prec_Trueness"] <-
            round(
                digits = 2,
                sum(na.rm = TRUE,
                    c(abs(median(a.df04[a.df04[, "Species"] == "Human"  , "log2FC"]) - expFC_human),
                      abs(median(a.df04[a.df04[, "Species"] == "Yeast"  , "log2FC"]) - expFC_yeast),
                      abs(median(a.df04[a.df04[, "Species"] == "E.coli"  , "log2FC"]) - expFC_ecoli),
                      abs(median(a.df04[a.df04[, "Species"] == "C.elegans"  , "log2FC"]) - expFC_celegans)
                    )))


        # Dispersion is the average distance of individual measurements around the respective median.
        summary_stats[, "Prec_Dispersion"] <-
            round(digits = 2,
                  mean(na.rm = TRUE,
                       c(
                           abs(a.df04[a.df04[, "Species"] == "Human"  , "log2FC"] - median(a.df04[a.df04[, "Species"] == "Human"  , "log2FC"])),
                           abs(a.df04[a.df04[, "Species"] == "Yeast"  , "log2FC"] - median(a.df04[a.df04[, "Species"] == "Yeast"  , "log2FC"])),
                           abs(a.df04[a.df04[, "Species"] == "E.coli"  , "log2FC"] - median(a.df04[a.df04[, "Species"] == "E.coli"  , "log2FC"])),
                           abs(a.df04[a.df04[, "Species"] == "C.elegans"  , "log2FC"] - median(a.df04[a.df04[, "Species"] == "C.elegans"  , "log2FC"]))
                       )))



        # Folder that contained the analysed pg and pc matrix
        summary_stats[, "folder_input"] <- folder_input




        return(summary_stats)

    }
