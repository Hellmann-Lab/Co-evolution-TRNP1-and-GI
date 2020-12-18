###runs as bash pipeline.sh <list of input files to process> <Phred quality to filter> <study name>
#loop through all data given as a list; weirdly not all arguments are part of that list ... echo -e allows escaping with backslash
echo -e "-------IGNORE-------"
#this counts all input to script that is contained in a list (here data files using a wildcard expression); it does curiously NOT include further parametes like Phred quality and output name ... they will produce an error message
records=$(ls $@ | wc -l)
echo -e "-------IGNORE-------\n"
echo -e "Pipeline applied to" $records "files.\n"
#above the length of a list (files to process) is stored in $records; in order for the loop only to run through these list elements we say "run through all input ($@) beginning with the 1st and ending with #$records" which is length of the list supplied
for files in ${@: 1:$records}
do
  #first remove directories (rm -r is recursive an also removes files within dictionaries) of files with the same name as the ones to be processed; else these directories can't be created ... still this should be kept in mind when naming sequencing files :P
  #rm -r 'filtered.'$files 
  #mkdir 'filtered.'$files
  #first (althouth strictly not necessary at all, because you only use it once) make shortcut for fqfilter function
  fqfilter=fqfilter.pl
  #count "reads" in files that are looped (looping variable $files), to count entries in FastQ count @ as every entry begins with @
  echo -e "Beginning...\n Counting records in original file\n"
  filterNone=$(zcat $files | grep -c @)
  #start filtering given file $files with given Phred value $(second to last value supplied to function) to the output specified in $(last value supplied to function); all tempory files are stored in the current working directory $PWD
  echo -e "starting the filtering process of" $filterNone "records\n with Phred >" ${@: -2:1} "\n output to" 'filtered.'$files/${@: -1}'.'$files'2.txt' "\n filtering report to" 'filtered.'$files/'report.'$files'2.txt'
  $fqfilter $files $files 1 ${@: -2:1} 1 ${@: -2:1} 1-5 6-10 12 tmp $PWD
  #when filtering is completed, display a message and give number of remaining data's lines; wc -l stores columns whitespace seperated, the second column holding information on the file processed, the first indicating #lines
  filterQual=$(wc -l tmp.barcodelist.filtered.sam | cut -d ' ' -f 1)
  #arithmetic operations are to be put within $((1+2)) or $(($a+$b))
  echo -e "\nquality filtering completed...\n remaining records" $filterQual "(removed" $(($filterNone-$filterQual)) "reads)"
  #remove unwanted fqfilter output; the names are not dynamic ... 
  rm tmp.cdnaread.filtered.fastq.gz tmp.barcoderead.filtered.fastq.gz
  #filter reads for constant gfp region; sequence information is stored in column 10 of sam files; the constant region length was chosen sort of arbitrary
  echo -e "\nnow filtering remaining reads against constant GFP region and cloning site"
  samtools view tmp.barcodelist.filtered.sam | cut -f 10 | grep -o -P '.{0,10}TCTAGAGTCGCGGCCTTACT' | sed 's/TCTAGAGTCGCGGCCTTACT//' | egrep '^.{10,10}$' > tmp.const.match2.txt
  #output a similar message to above
  filterConst=$(wc -l tmp.const.match2.txt | cut -d ' ' -f 1)
  echo -e "\nconstant region filtering completed...\n remaining records" $filterConst "(removed" $(($filterQual-$filterConst)) "reads)"
  #extract only barcodes and count their occurences (in order to not store huge files in working memory, uniq only looks at adjacent lines, so data needs to be sorted first - sort -n -r sorts numerically, beginning with the largest); sent the output to the specified file
  cut -c 1-10 tmp.const.match2.txt | sort | uniq -c | sort -n -r > 'filtered.'$files/${@: -1}'.'$files'2.txt'
  #remove temporary data
  #rm tmp.barcodelist.filtered.sam tmp.const.match.txt
  #display data preview; barcodes are now reverse complement of the information stored in the original library file
  echo -e "\n-----it's done, this is what it looks like-----\n"
  head 'filtered.'$files/${@: -1}'.'$files'2.txt'
  #store read development into an out file where > creates this file with the first input an >> appends all following inputs to the created file
  #the way this is written, a whitespace is introduced in front of the "applied filter"
  echo -e $filterNone '\t' "total reads" > 'filtered.'$files/'report.'$files'2.txt'
  echo -e $filterQual '\t' "quality filtered reads"  >> 'filtered.'$files/'report.'$files'2.txt'
  echo -e $filterConst '\t' "constant region filtered reads"  >> 'filtered.'$files/'report.'$files'2.txt'
  echo -e "\n"
  cat 'filtered.'$files/'report.'$files'2.txt'
done
