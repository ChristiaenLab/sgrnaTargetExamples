#get targets
java -Xmx4g -jar FlashFry-assembly-1.12.jar \
 discover \
 --database kh2013_db \
 --fasta knowntargets.fasta \
 --positionOutput \
 --output knowntargets.output

#get scores
java -Xmx4g -jar FlashFry-assembly-1.12.jar \
 score \
 --input knowntargets.output \
 --output knowntargets.output.scored \
 --scoringMetrics doench2014ontarget,doench2016cfd,dangerous,hsu2013,minot \
 --database kh2013_db \

bedtools getfasta -name -fi JoinedScaffold.fasta -fo targetSites.fasta -bed targetSites.bed -s

#get targets
java -Xmx4g -jar FlashFry-assembly-1.12.jar \
 discover \
 --database kh2013_db \
 --fasta targets.fasta \
 --positionOutput \
 --output targets.output

#get scores
java -Xmx4g -jar FlashFry-assembly-1.12.jar \
 score \
 --input targets.output \
 --output targets.output.scored \
 --scoringMetrics doench2014ontarget,doench2016cfd,dangerous,hsu2013,minot \
 --database kh2013_db \

