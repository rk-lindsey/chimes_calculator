#!/bin/bash
for File in *.0.2b_clusters.txt
do
Name=$(basename $File)
Frame=$(echo $Name | cut -d '.' -f1)
cat ${Frame}.*.2b_clusters.txt > ${Frame}.all-2b-clusters.txt
rm ${Frame}.*.2b_clusters.txt
rm ${Frame}.*.2b_list.txt
done
for File in *.0.3b_clusters.txt
do
Name=$(basename $File)
Frame=$(echo $Name | cut -d '.' -f1)
cat ${Frame}.*.3b_clusters.txt > ${Frame}.all-3b-clusters.txt
rm ${Frame}.*.3b_clusters.txt
rm ${Frame}.*.3b_list.txt
done
for File in *.0.4b_clusters.txt
do
Name=$(basename $File)
Frame=$(echo $Name | cut -d '.' -f1)
cat ${Frame}.*.4b_clusters.txt > ${Frame}.all-4b-clusters.txt
rm ${Frame}.*.4b_clusters.txt
rm ${Frame}.*.4b_list.txt
done
