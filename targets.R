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

snps <- import("SNPs.gff3")
findOverlaps(coords,snps)

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

overlapsAny(coords,targetSites)

export(targetSites,'targets.bed')

system2('bash','flashfry.sh')

#split matches for each sequence
parseOutput <- function(x) {
	seqnames <- sub('.*<(.*):.*','\\1',x)
	strand <- sapply(grepl('\\^R',x), function(x) if(x) '-' else '+')
	coord <- as.numeric(sub('.*:([0-9]+).*','\\1',x))
	ranges <- IRanges(coord+1,coord+23)
	mismatches <- sub('.*_[0-9]+_([0-9]+)<.*','\\1',x)
	sequence <- sub('_.*','',x)
	return(GRanges(
		       seqnames,ranges,strand,
		       mismatches=as.numeric(mismatches),
		       sequence=sequence))
}

#select 3 best targets for each gene
best3 <- function(x) {
	x <- x[order(x$SNPct)]
	x <- x[order(x$otCount)]
	return(x[1:min(length(x),3)])
}

out <- read.delim('targets.output')
outscore <- read.delim('targets.output.scored')
out <- merge(out,outscore)
out <- out[out$Doench2014OnTarget>0.65&is.finite(out$Doench2014OnTarget),]

offTargets <- strsplit(out$offTargets,',')

offTargets <- do.call(GRangesList,lapply(offTargets, parseOutput))

#select which match is on-target for each seq
onTarget <- unlist(
	do.call(
		GRangesList,
		sapply(
		       offTargets,
		       function(x) x[order(x$mismatches)[1]]
		)
	)
)

mcols(onTarget) <- cbind(out[,-2:-4],mcols(onTarget))

onTarget <- onTarget[!duplicated(onTarget)]

#look for SNPs in each on-target
onTarget$SNPct <- countOverlaps(onTarget,snp)

#split seqs by gene
onTarget <- split(onTarget,sub('\\.v.*','',onTarget$contig))

#select 3 best seqs for each gene
onTarget <- lapply(onTarget,best3)
onTarget <- unlist(do.call(GRangesList,onTarget))

write.csv(mcols(onTarget),'targetSiteScore.csv')


