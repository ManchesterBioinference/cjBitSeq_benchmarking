unzip genes.gtf.zip
cd WholeGenomeFasta
unzip genome.fa.zip
bowtie2-build -f genome.fa genome
cd ..
tophat2 -G genes.gtf --transcriptome-index=transcriptome_data/known WholeGenomeFasta/genome && rm -r tophat_out
grep -e '^2L\|^2R\|^3L\|^3R\|^4\|^M\|^X' genes.gtf > genesNew.gtf
grep -v 'Het' genesNew.gtf > tmp.gtf
mv tmp.gtf genesNew.gtf
mkdir rsemReference
cd rsemReference
rsem-prepare-reference ../transcriptome_data/known.fa ref
bowtie2-build -f ref.transcripts.fa ref

