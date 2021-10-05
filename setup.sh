#download KH genome sequencf
wget -U firefox http://ghost.zool.kyoto-u.ac.jp/datas/JoinedScaffold.zip
unzip -o JoinedScaffold.zip

#download SNPs
wget -U firefox http://ghost.zool.kyoto-u.ac.jp/cgi-bin/gb2/gbrowse/kh/?l=SNPS;f=save+datafile

#download flashfry
wget https://github.com/mckennalab/FlashFry/releases/download/1.12/FlashFry-assembly-1.12.jar

#create db
mkdir -p build
java -Xmx4g -jar flashfry/FlashFry-assembly-1.12.jar \
 index \
 --tmpLocation ./build \
 --database kh2013_db \
 --reference JoinedScaffold \
 --enzyme spcas9ngg
