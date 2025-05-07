library(tidyverse)
library(SingleCellExperiment)
library(edgeR)
library(batchelor)
library(wrMisc)
library(pheatmap)
library(RColorBrewer)

# Panel 1: Expression of lymphoid, DAM and homeostatic genes in human homeostatic and DAM clusters
# Heatmap colors represent row-scaled, normalized and batch corrected mean expression values per cluster. 

# make pseudobulks per cluster
# the input sce is a SingleCellExperiment object created from the Green et al data (syn31512863), subset to myeloid cells 

make_cluster_pseudobulks <- function(sce) { 

    ## remove samples from -B6 and -alone batches, most samples run in duplicate "A" and "B" batches
    f<-grepl(sce$library_batch, pattern="-B6$")|grepl(sce$library_batch,pattern="alone$")
    sce <- sce[,!f] 

    ### fix erroneous batch info so that prefix is the same for A and B batches 
    f<-sce$library_batch=="201024-B59-B"
    sce$library_batch[f] <- "201014-B59-B"

    ### combine counts across A and B batches for each individual
    ### first set the grouping info
    batch <-strsplit(sce$library_batch,"-") %>% sapply(function(x) x[1])
    id <- strsplit(sce$library_batch,"-") %>% sapply(function(x) x[2])
    a_vs_b <- strsplit(sce$library_batch,"-") %>% sapply(function(x) x[3])
    sce$new_batch <- paste(batch,id,sep="-")
    sce$celltype_donor_newbatch <- paste(sce$state, sce$individualID, sce$new_batch, sep="_")

    ### next create pseudobulk and coldata
    pseudo_bulk <- t(rowsum(t(counts(sce)), sce$celltype_donor_newbatch))
    
    cdat <- as.data.frame(colData(sce)) %>% 
    group_by(celltype_donor_newbatch) %>% 
    mutate(ncells=n()) %>%
    distinct() %>% 
    ungroup()
    
    cdat <- cdat[match(colnames(pseudo_bulk), cdat$celltype_donor_newbatch),]

    ### turn it into sce
    pb = SingleCellExperiment(assays=list(counts=pseudo_bulk),
                              rowData = rowData(sce),
                              colData = as(cdat,"DataFrame"))

    rownames(pb) <- rowData(pb)$symbol 

    ### remove pseudobulks with <10 cells 
    pb <- pb[,pb$ncells>=10]

    return(pb)
}

pb <- make_cluster_pseudobulks(sce) ## 2057 total samples

# remove duplicates per donor -- some run in multiple batches, not just A/B

    keep <- as.data.frame(colData(pb)) %>% group_by(donor_state) %>% slice_max(order_by=ncells) %>% pull(celltype_donor_newbatch)
    f<- pb$celltype_donor_newbatch %in% keep ## 2017, still one tie
    pb <- pb[,f]
    f <- pb$celltype_donor_newbatch=="Mic.3_R8760165_200825-B48"
    pb <- pb[,!f] ## now 2016 samples 

# remove macrophages and monocytes

    f<- pb$cell.type %in% c("Macrophages","Monocytes")
    pb <- pb[,!f] ## 1963

# batch correct for plotting 

    pb2 <- edgeR::calcNormFactors(pb)
    pb2 <- edgeR::cpm(pb2, log=TRUE)
    pb2b <- batchelor::regressBatches(pb2, batch=pb$new_batch)
    rownames(pb2b) <- rownames(pb)
    assays(pb)$corrected <- assays(pb2b)$corrected
    rownames(pb) <- rowData(pb)$symbol
    rownames(pb2) <- rownames(pb)
    assays(pb)$logcounts <- pb2

# add annotations 

    pb$cogdx_level <- ifelse(pb$cogdx<2,"nci",ifelse(pb$cogdx<4,"mci",ifelse(pb$cogdx<6,"ad", "other")))
    pb$cogdx_level <- factor(pb$cogdx_level, levels=c("nci","mci","ad"))

# plot 

    lymphoid <- c("CD28","PDCD1","IL2RG","DKK2","CD6","CD72","CXCR4")
    dam <- c("APOE","CD163","SPP1","CD9","LPL")
    homeo <- c("CSF1R","TMEM119","CX3CR1","CD33","P2RY12")
    gs <- c(homeo, lymphoid, dam)

    annotation_df <- data.frame(gene=gs, group=NA) 
    f<- annotation_df$gene %in% homeo
    annotation_df$group[f] <- "homeo"
    f<- annotation_df$gene %in% dam
    annotation_df$group[f] <- "dam"
    f<- annotation_df$gene %in% lymphoid
    annotation_df$group[f] <- "lymphoid"
    rownames(annotation_df) <- gs

    annot_column <-  as.data.frame(colData(pb)[,c("state", "state_group")])

    mat_pb <- wrMisc::rowGrpMeans(as.matrix(assays(pb)$logcounts), pb$state)
    rownames(mat_pb) <- rowData(pb)$symbol
    annot_col <- annot_column %>% select(state, state_group) %>% distinct()
    rownames(annot_col) <- annot_col$state
    f<- annotation_df$group %in% c("homeo","lymphoid","dam")
    gs <- annotation_df$gene[f]

    filt <- annot_col$state_group %in% c("lipid_assoc","surveilling") & annot_col$state!="Mic.5"

    pheatmap::pheatmap(mat_pb[gs,filt], cluster_rows=TRUE, cluster_cols=TRUE,
    color = colorRampPalette(rev(brewer.pal(n = 7, name =
    "RdBu")))(100),scale="row",breaks=seq(-2,2,length=101),
    annotation_col=annot_col[,c("state_group"), drop=F],
    annotation_row=annotation_df[f,"group", drop=F],
    fontsize=8, show_rownames=TRUE, show_colnames=TRUE, border_color="white")

# Panel 2: Expression of SPI1 in PU.1-risk carriers with and without AD
# Normalized SPI1 expression values are plotted for each sample according to their cognitive diagnosis and PU.1 risk variant status (rs1057233 A/A = risk, G/A and G/G = non-risk). No significant differences are observed between AD and control, or between risk and non-risk carriers, although risk carriers have slightly higher SPI1 levels in AD consistent with eQTL results in monocytes.

# make pseudobulks per donor
# the input sce is a SingleCellExperiment object created from the Green et al data (syn31512863)

    make_pseudobulk <- function(sce) { 
    ## get data, remove samples from -B6 and -alone batches
    f<-grepl(sce$library_batch, pattern="-B6$")|grepl(sce$library_batch,pattern="alone$")
    sce <- sce[,!f] 

    ### fix B59 so that prefix is the same for A and B batches 
    f<-sce$library_batch=="201024-B59-B"
    sce$library_batch[f] <- "201014-B59-B"

    ### combine counts across A and B batches for each individual
    ### first set the grouping info
    batch <-strsplit(sce$library_batch,"-") %>% sapply(function(x) x[1])
    id <- strsplit(sce$library_batch,"-") %>% sapply(function(x) x[2])
    a_vs_b <- strsplit(sce$library_batch,"-") %>% sapply(function(x) x[3])
    sce$new_batch <- paste(batch,id,sep="-")
    sce$celltype_donor_newbatch <- paste(sce$cell.type, sce$individualID, sce$new_batch, sep="_")

    ### next create pseudobulk and coldata
    pseudo_bulk <- t(rowsum(t(counts(sce)), sce$celltype_donor_newbatch))
    
    cdat <- as.data.frame(colData(sce)) %>% 
    group_by(celltype_donor_newbatch) %>% 
    mutate(ncells=n()) %>%
    distinct() %>% 
    ungroup()
    
    cdat <- cdat[match(colnames(pseudo_bulk), cdat$celltype_donor_newbatch),]

    ### turn it into sce
    pb = SingleCellExperiment(assays=list(counts=pseudo_bulk),
                              rowData = rowData(sce),
                              colData = as(cdat,"DataFrame"))

    rownames(pb) <- rowData(pb)$symbol 

    ### remove pseudobulks with <10 cells 
    pb <- pb[,pb$ncells>=10]

    return(pb)
    }

    pb <- make_pseudobulk_rosmap(sce)

    ### keep microglia only 
    f<- pb$cell.type=="Microglia"
    pb <- pb[,f] ## 440 

    ### remove duplicates per donor -- some run in multiple batches, not just A/B
    keep <- as.data.frame(colData(pb)) %>% group_by(individualID) %>% slice_max(order_by=ncells) %>% pull(celltype_donor_newbatch)
    f<- pb$celltype_donor_newbatch %in% keep ## 428 
    pb <- pb[,f]

# add donor genotypes - need access to genotypes from synapse, match on projid for rs1057233
# recode genotypes to risk/non-risk
    sce$PU1_risk <- ifelse(sce$genotype == "A/A", "risk", "non-risk")

# normalize 
    pb2 <- edgeR::calcNormFactors(pb)
    pb2 <- edgeR::cpm(pb2, log=TRUE)
    assays(pb)$logcounts <- pb2
    rownames(pb) <- rowData(pb)$symbol
    rownames(pb2) <- rownames(pb)

# add annotations & SPI1 logCPM to coldata
    pb$cogdx_level <- ifelse(pb$cogdx<2,"nci",ifelse(pb$cogdx<4,"mci",ifelse(pb$cogdx<6,"ad", "other")))
    pb$cogdx_level <- factor(pb$cogdx_level, levels=c("nci","mci","ad"))
    pb$SPI1_lc <- assays(pb)$logcounts["SPI1",]

# create dataframe for plotting

    to_plot <- as.data.frame(colData(pb)) %>% select(cogdx_level, SPI1_lc, PU1_risk)

# plot 

    to_plot %>% 
    ggplot(aes(x=cogdx_level,y=SPI1_lc, col=PU1_risk))+stat_summary()+theme_classic()+
    labs(x="Cognitive diagnosis",y="SPI1 logCPM")+ggtitle("SPI1 expression vs. PU.1 carrier status")+
    scale_color_manual(values=c("black","tomato"))+labs(col="PU1 carrier\nstatus")+
    theme(axis.text.x=element_text(color="black",size=11))+theme(axis.text.y=element_text(color="black",size=11))

# Panel 3: PU.1 risk carriers have fewer DAM in AD 
# the input sce is a SingleCellExperiment object created from the Green et al data (syn31512863) with PU1 risk status and cogdx_level added as shown above

# calculate state % per donor

    f <- sce$cell.type=="Microglia"
    sce <- sce[,f]

    x <- as.data.frame(colData(sce)) %>% group_by(projid) %>% 
                        summarize(total=n(),
                        mic13_all=sum(state=="Mic.13")/n())

# create dataframe for plotting

    to_plot <- left_join(x, as.data.frame(colData(sce)) %>% select(cogdx_level, mic13_all, PU1_risk))

# plot 

    to_plot %>% 
    ggplot(aes(x=cogdx_level,y=mic13_all, col=PU1_risk))+stat_summary()+theme_classic()+
    labs(x="Cognitive Diagnosis",y="Proportion of DAM microglia (Mic13)")+
    ggtitle("Mic13 cellularity")+labs(col="PU.1 carrier\nstatus")+
    scale_color_manual(values=c("black","tomato"))+
    theme(axis.text.x=element_text(color="black",size=11))+
    theme(axis.text.y=element_text(color="black",size=11))
    