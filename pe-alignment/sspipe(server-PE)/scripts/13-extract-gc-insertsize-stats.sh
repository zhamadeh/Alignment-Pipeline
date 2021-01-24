#Calculates mode GC content and insert size for each library. Some name processing will be sample name specific (i.e. the sed statements and the paste line may need to be edited)
#Run like this: bash extract_gc_insertsize_stats.sh /path/to/insert/size/directory /path/to/GC/metrics/directory /path/to/output

#############################################################################################################################################################################################

#Make a list of sample names
ls $1/*.colinsert.txt | sed 's/\.trimmed.colinsert\.txt//g' | rev | cut -f1 -d/ | rev > list


#Initiate an ouput file and print headers
printf "Library Mode_GC Mode_insert_size\n" > tmp1


#Loop over sample names and extract mode GC and mode insert size, then print to file
while read name;
do

        paste <(echo $name | awk -F / '{print $NF}' ) <(grep -v '^#' $2/$name.trimmed.gc_bias_metrics.txt | awk '{print $4,$6}' | sort -k2,2nr | head -1 | cut -f1 -d" ") <(grep -A1 'MODE_INSERT_SIZE' $1/$name.trimmed.colinsert.txt | tail -1 | cut -f2) >> tmp1


done <list

awk '{print $1,$2,$3}' OFS="\t" tmp1 > $3/gc_insertsize_stats.txt

rm list tmp1
