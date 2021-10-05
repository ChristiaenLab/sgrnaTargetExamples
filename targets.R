library(GenomicRanges)
library(rtracklayer)
library(Biostrings)
library(BSgenome.Cintestinalis.KH.JoinedScaffold)
library(GenomicFeatures)

dat <- read.csv('primers.csv')

targets <- DNAStringSet(gsub('\\W','',dat$Protospacer...PAM.Sequence..protospacer..specific..20bp.before.PAM.in.genome.))

coords <- dat$Coordinates.of.the.sgRNA
coords <- sub('\\.\\.','\\-',gsub(',','',coords))

names(targets) <- paste0(dat$sgRNA.name,'::',coords)
export(targets,'knowntargets.fasta')

coords <- do.call(rbind,strsplit(coords,':|\\-'))
coords <- GRanges(
	   coords[,1],
	   IRanges(as.numeric(coords[,2]),as.numeric(coords[,3]))
)
names(coords) <- dat$sgRNA.name

system2('./flashfry.sh')

kh <- makeTxDbFromGFF('KH.KHGene.2013.gff3')

tpt <- transcripts(kh)
promoter <- promoters(tpt,1107,107)
exon <- exonsBy(kh,'tx')
names(exon) <- tpt$tx_name

#check order for exons and transcripts
all(unlist(unique(seqnames(exon)))==seqnames(tpt))

intron <- intronsByTranscript(kh)
names(intron) <- tpt$tx_name

promotersel <- unlist(sapply(paste0(unique(dat$KH.gene.ID),'.v'),grep,promoter$tx_name))
targetSite <- promoter[promotersel]
names(targetSite) <- paste0(targetSite$tx_name,'_promoter')

exonsel <- unlist(sapply(paste0(unique(dat$KH.gene.ID),'.v'),grep,names(exon)))
targetExons <- exon[exonsel]
targetExons <- unlist(do.call(GRangesList,sapply(targetExons,function(x) x[1:min(3,length(x))])))
names(targetExons) <- paste0(names(targetExons),'_exon',as.character(targetExons$exon_rank))

targetSites <- c(targetSite,targetExons)

tmp <- reduce(targetSites)
findOverlaps(tmp,tpt)

export(targetSites,'flashfry/targets.bed')


