#!/bin/bash
#for str in "FBtr0073806"
#do
#echo $str
for i in *.cufflinks.out
do
cd $i
cp transcripts.gtf transcriptsCopy.gtf
awk '!/FBtr0079998/' transcriptsCopy.gtf > temp && mv temp transcriptsCopy.gtf
awk '!/FBtr0079997/' transcriptsCopy.gtf > temp && mv temp transcriptsCopy.gtf
awk '!/FBtr0073570/' transcriptsCopy.gtf > temp && mv temp transcriptsCopy.gtf
awk '!/FBtr0073806/' transcriptsCopy.gtf > temp && mv temp transcriptsCopy.gtf
cd ../
echo $i
done
#done

