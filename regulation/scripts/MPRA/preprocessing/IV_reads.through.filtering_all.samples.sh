#borrows most functionality from script concatenating fastqc information

rm ../read_summaries/reads.through.filtering_all.samples2.txt
touch ../read_summaries/reads.through.filtering_all.samples2.txt
echo -e "reads\tstep\tsample\treplicate" > ../read_summaries/reads.through.filtering_all.samples2.txt #add Sample column header to header line and header line to final output file also remove # from #Base so it can be read by R
for file in $@ #loop through input files
do
  sample=$(echo $file | cut -d '.' -f2) #extract sample from file name, sample needs to be seperated by _
replicate=$(echo $file | cut -d '-' -f3,4) #extract sample from file name, sample needs to be seperated by _
  cat $file | while read  line #pipe temp file data into while loop, read each line and store it in <line> to be modified; here for all table contents add the correct sample name in the final column seperated by a \t 
    do echo -e "$line\t$sample\t$replicate" >> ../read_summaries/reads.through.filtering_all.samples2.txt
    done
done