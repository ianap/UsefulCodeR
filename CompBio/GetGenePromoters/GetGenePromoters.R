library('GenomicRanges')
library('GenomicFeatures')
library('TxDb.Hsapiens.UCSC.hg19.knownGene')

#biocLite("ChIPpeakAnno")
library('ChIPpeakAnno')
library(biomaRt)
library(org.Hs.eg.db)
library(data.table)
library(compare)

#read annotation file from BETA hg19
refseq_hg19<-read.table("hg19.refseq")
refseq_hg19<-as.character(refseq_hg19[,1])

#mapping entrezid to refseq
my_keys <- refseq_hg19
a <- select(org.Hs.eg.db,
            keys = my_keys,
            columns=c("REFSEQ","ENTREZID", "SYMBOL"),
            keytype="REFSEQ")

b <- unique(select(TxDb.Hsapiens.UCSC.hg19.knownGene,
            keys = a$ENTREZID,
            columns=c('GENEID', 'TXCHROM', 'TXSTART', 'TXEND', 'TXID'),
            keytype="GENEID"))

names(b) <- c('ENTREZID', 'TXID', 'TXCHROM', 'TXSTART', 'TXEND')
c <- (unique(merge(a, b, 'ENTREZID'))
c<-c[complete.cases(c),]

# Genes and exons coordinates
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

genes <- as.data.table(as.data.frame(transcripts(txdb)))
exs <- as.data.table(as.data.frame(exons(txdb)))

#set TSS, if strand "+" then TSS=startGene, else TSS = endGene
genes[,tss:=ifelse(strand=="+",start,end) ]
promoters<-as.data.table(as.data.frame(promoters(txdb)))

genes_pos_strand<-subset(genes,strand=="+")
genes_neg_strand<-subset(genes,strand=="-")

promoters_positive_strands<-subset(promoters, strand=="+")
promoters_negative_strands<-subset(promoters, strand=="-")

col_names<-c("chr","start","end")
promoters_positive_strands<-unique(promoters_positive_strands[,1:3, with=FALSE])
promoters_negative_strands<-unique(promoters_negative_strands[,1:3, with=FALSE])
exons<-unique(exs[,1:3,with=FALSE])

colnames(promoters_positive_strands)<-col_names
colnames(promoters_negative_strands)<-col_names
colnames(exons)<-col_names

##### For a particular transcription factor get the    #####
##### number of TF's peaks up- downstream the promoter.#####

#read ChiP-seq
tf<-data.table(read.table("IRF4.txt",col.names = col_names))
#consider peak center only 
tf[,pcStart:=as.integer((start+end)/2)]
tf[,pcEnd:=pcStart+1]
#remove peak  start and peak end
tf<-tf[,c(1,4,5), with=FALSE]
colnames(tf)<-c("chr","start","end")
setkey(tf, chr, start, end)

overlap_positive<-foverlaps(promoters_positive_strands, tf, type="any", which=TRUE)
overlap_positive<-overlap_positive[complete.cases(overlap_positive),]
peaks_in_pos_promoters<-unique(tf[as.numeric(overlap_positive$yid),])

overlap_negative<-foverlaps(promoters_negative_strands,tf,type="any",which=TRUE)
overlap_negative<-overlap_negative[complete.cases(overlap_negative)]
peaks_in_neg_promoters<-unique(tf[as.numeric(overlap_negative$yid),])

overlap_exons<-foverlaps(exons, tf, type="any", which=TRUE)
overlap_exons<-overlap_exons[complete.cases(overlap_exons),]
peaks_in_exons<-unique(tf[as.numeric(overlap_exons$yid),])

#getting up- and downstream regions
#1. Downn- and Upstream positive strand
promoter_pos_strand_tss<-unique(promoters_positive_strands)
promoter_pos_strand_tss[,tss:=(start+2000) ]

promoter_pos_strand_down<-promoter_pos_strand_tss[,c(1,3,4),with=FALSE]
setcolorder(promoter_pos_strand_down, c("chr", "tss", "end"))
colnames(promoter_pos_strand_down)<-col_names
ovlap_pos_prom_down<-foverlaps(promoter_pos_strand_down,tf,type="any",which=TRUE)
ovlap_pos_prom_down<-ovlap_pos_prom_down[complete.cases(ovlap_pos_prom_down),]
peaks_in_pos_prom_down<-unique(tf[as.numeric(ovlap_pos_prom_down$yid),])

promoter_pos_strand_up<-promoter_pos_strand_tss[,c(1,2,4),with=FALSE]
setcolorder(promoter_pos_strand_up, c("chr", "start", "tss"))
colnames(promoter_pos_strand_up)<-col_names
ovlap_pos_prom_up<-foverlaps(promoter_pos_strand_up,tf,type="any",which=TRUE)
ovlap_pos_prom_up<-ovlap_pos_prom_up[complete.cases(ovlap_pos_prom_up),]
peaks_in_pos_prom_up<-unique(tf[as.numeric(ovlap_pos_prom_up$yid),])

#3. Upstream negative strand
promoter_neg_strand_tss<-unique(promoters_negative_strands)
promoter_neg_strand_tss[,tss:=(end-2000) ]

promoter_neg_strand_down<-promoter_neg_strand_tss[,c(1,2,4),with=FALSE]
setcolorder(promoter_neg_strand_down, c("chr", "start", "tss"))
colnames(promoter_neg_strand_down)<-col_names
ovlap_neg_prom_down<-foverlaps(promoter_neg_strand_down,tf,type="any",which=TRUE)
ovlap_neg_prom_down<-ovlap_neg_prom_down[complete.cases(ovlap_neg_prom_down),]
peaks_in_neg_prom_down<-unique(tf[as.numeric(ovlap_neg_prom_down$yid),])

promoter_neg_strand_up<-promoter_neg_strand_tss[,c(1,3,4),with=FALSE]
setcolorder(promoter_neg_strand_up, c("chr", "tss", "end"))
colnames(promoter_neg_strand_up)<-col_names
ovlap_neg_prom_up<-foverlaps(promoter_neg_strand_up,tf,type="any",which=TRUE)
ovlap_neg_prom_up<-ovlap_neg_prom_up[complete.cases(ovlap_neg_prom_up),]
peaks_in_neg_prom_up<-unique(tf[as.numeric(ovlap_neg_prom_up$yid),])
#overlap peaks and promoers up and down:
setkey(tf, chr, start, end)

overlap_positive_up<-foverlaps(promoter_pos_strand_up, tf, type="any", which=TRUE)
overlap_positive_up<-overlap_positive_up[complete.cases(overlap_positive),]
peaks_in_pos_promoters<-unique(tf[as.numeric(overlap_positive_up$yid),])

#compare two dataframes
require(sqldf)
intersection <- sqldf('SELECT * FROM peaks_in_pos_promoters INTERSECT SELECT * FROM peaks_in_exons')
intersection <- sqldf('SELECT * FROM peaks_in_pos_promoters INTERSECT SELECT * FROM peaks_in_neg_promoters')
