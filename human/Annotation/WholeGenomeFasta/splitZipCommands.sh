split --bytes=150M genome.fa genomeSplitted

for s in genomeSplitted*; do zip zipped${s}.zip ${s}; done

rm genome.fa 
