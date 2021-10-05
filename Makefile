FASTA = JoinedScaffold
URL = http://ghost.zool.kyoto-u.ac.jp/datas/$(FASTA).zip
DB = kh2013_db
SNPS = http://ghost.zool.kyoto-u.ac.jp/cgi-bin/gb2/gbrowse/kh/?l=SNPS;f=save+datafile
GENOME = KH.KHGene.2013.gff3
VERSION = 1.12

install: build
	java -Xmx4g -jar FlashFry/FlashFry-assembly-$(VERSION).jar \
	 index \
	 --tmpLocation ./build \
	 --database $(DB) \
	 --reference $(FASTA) \
	 --enzyme spcas9ngg

build: FlashFry-assembly-$(VERSION).jar
	mkdir -p build

FlashFry-assembly-$(VERSION).jar:
	wget https://github.com/mckennalab/FlashFry/releases/download/$(VERSION)/FlashFry-assembly-$(VERSION).jar

$(GENOME):
	wget -U firefox http://ghost.zool.kyoto-u.ac.jp/datas/$(GENOME).zip
	unzip $(GENOME).zip

SNPS.gff3:
	wget -U firefox $(SNPS)

$(FASTA):
	wget -U firefox $(URL)
	unzip -o $(FASTA).zip

clean:
	rm -rf build
	rm -f $(GENOME).zip
	rm -f $(FASTA).zip
