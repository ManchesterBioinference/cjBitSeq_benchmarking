unzip genes.gtf.zip
cd WholeGenomeFasta
for s in zipped*; do unzip ${s}; done
cat genomeSplitted* > genome.fa
bowtie2-build -f genome.fa genome
cd ..
tophat2 -G genes.gtf --transcriptome-index=transcriptome_data/known WholeGenomeFasta/genome && rm -r tophat_out

grep -e '^chrM\|^chr1\|^chr2\|^chr3\|^chr4\|^chr5\|^chr6\|^chr7\|^chr8\|^chr9\|^chr10\|^chr11\|^chr12\|^chr13\|^chr14\|^chr15\|^chr16\|^chr17\|^chr18\|^chr19\|^chr21\|^chr22\|^chrX\|^chrY' genes.gtf > genesNew.gtf
# when you run the previous command, is not able to remove lines with chr17_ctg5_hap1. So type next
grep -v 'chr17_ctg5_hap1' genesNew.gtf > tmp.gtf
mv tmp.gtf genesNew.gtf
# when you run the previous command, is not able to remove lines with chr17_gl000205_random. So type next
grep -v 'chr17_gl000205_random' genesNew.gtf > tmp.gtf
mv tmp.gtf genesNew.gtf
grep -v 'chr19_gl000209_random' genesNew.gtf > tmp.gtf
mv tmp.gtf genesNew.gtf
grep -v 'chr1_gl000191_random' genesNew.gtf > tmp.gtf
mv tmp.gtf genesNew.gtf
grep -v 'chr1_gl000192_random' genesNew.gtf > tmp.gtf
mv tmp.gtf genesNew.gtf
grep -v 'chr4_ctg9_hap1' genesNew.gtf > tmp.gtf
mv tmp.gtf genesNew.gtf
grep -v 'chr4_gl000193_random' genesNew.gtf > tmp.gtf
mv tmp.gtf genesNew.gtf
grep -v 'chr4_gl000194_random' genesNew.gtf > tmp.gtf
mv tmp.gtf genesNew.gtf
grep -Ev 'hap1|random' genesNew.gtf > tmp.gtf
mv tmp.gtf genesNew.gtf
grep -Ev 'hap' genesNew.gtf > tmp.gtf
mv tmp.gtf genesNew.gtf



mkdir rsemReference
cd rsemReference
rsem-prepare-reference ../transcriptome_data/known.fa ref
bowtie2-build -f ref.transcripts.fa ref







