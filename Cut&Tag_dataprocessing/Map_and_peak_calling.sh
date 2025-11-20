sn=$1 # sample name
$out_dir="raw_data"
R1="unspecified/*.UnspecifiedIndex.*.s_1.r_1.fq.gz"
R2="unspecified/*.UnspecifiedIndex.*.s_1.r_2.fq.gz"

# Step 1: Demultiplex the sequencing data
demuxFQ -b $1.unmatched_barcodes.r1.fq -c -d -i -e -t 1 -r 0.01 -R -l 9 -o $out_dir -s demul_summary_r1.txt demul_index_r1 $R1
demuxFQ -b $1.unmatched_barcodes.r2.fq -c -d -i -e -t 1 -r 0.01 -R -l 9 -o $out_dir -s demul_summary_r2.txt demul_index_r2 $R2




# Step 2: quality control
fastqc -o fastqc_results raw_data/${sn}r1.fq.gz raw_data/${sn}r2.fq.gz
multiqc fastqc_results #view the QC file




# Step 3: trim the sequencing data 
for file in *.r1.*.fastq
do
fq1=$file
fq2=${fq1/r1/r2}
ls $fq1
ls $fq2
echo "----"
sbatch --time 01:00:00 -o %j.out -e %j.err --mem 16000 --wrap "cutadapt -q 20  -o trimmed/${fq1%%.fastq}.trimmed.fq.gz -p trimmed/${fq2%%.fastq}.trimmed.fq.gz $fq1 $fq2 "
done




# Step 4: quality control
fastqc -o trimmed_fastqc_results trimmed/${sn}r1.trimmed.fq.gz raw_data/${sn}r2.trimmed.fq.gz
multiqc trimmed_fastqc_fastqc_results #view the QC file





# Step5: alignment with bwa
mkdir aligned
g='path/mm10_selected_ecoli_K12.fa'
w='path/mm10_eco.whitelist.sorted.bed'
  
#aln with bwa
for file in *.r1*trimmed.fq.gz 
  do
  f1=$file
  f2=${file/r1/r2}
  sbatch --time 12:00:00  -p epyc --mem 16G --wrap "bwa mem -M $g $f1 $f2  | samtools view -Sb -F780 -q 10 -L $w - | samtools sort -@ 8 -   > aligned/${f1%%.trimmed.fq.gz}.mm10.sort.bam"
done




# Step 6: Merge
for file in *.r1.00.mm10.sort.bam
do
sample=${file%%.r1.00.mm10.sort.bam}
sbatch --time 01:00:00 -p epyc --mem 8G --wrap "samtools merge $sample.merged.bam ${file%%.r1.00.mm10.sort.bam}.r1.*.mm10.sort.bam"
echo "-------"
done



# Step 7: mark duplicates with Picard
for bam in *.bam
do
sbatch -p epyc --time 12:00:00 --mem 8G --wrap "java -Xmx7g -jar /home/simeon01/applications/picard-2.20.3.jar MarkDuplicates INPUT=$bam OUTPUT=${bam%%.bam}.markduplicates.bam REMOVE_DUPLICATES=true  AS=true METRICS_FILE=${bam%%.bam}.markduplicates_metrics.txt | samtools index ${bam%%.bam}.markduplicates.bam"
done






# Step 8: generate count per million reads (CPM) normalized profile #play with diff fragment sizes
for file in *merged*.bam
do
sbatch --time 02:00:00 --mem 4G  -p epyc --wrap "bamCoverage --bam $file -o ${file%%bam}.150.cpm.bs10.bw --binSize 10 --normalizeUsing CPM  --extendReads --maxFragmentLength 150"
sbatch --time 02:00:00 --mem 4G  -p epyc --wrap "bamCoverage --bam $file -o ${file%%bam}.120.cpm.bs10.bw --binSize 10 --normalizeUsing CPM  --extendReads --maxFragmentLength 120"
sbatch --time 02:00:00 --mem 4G  -p epyc --wrap "bamCoverage --bam $file -o ${file%%bam}.130.cpm.bs10.bw --binSize 10  --normalizeUsing CPM  --extendReads --maxFragmentLength 130"
done








# For Peak callinf (stringent and relaxed mode)
for file in *merged*.bam
do

sbatch --time 12:00:00 --mem 4G   --wrap "bamCoverage --bam $file -o ${file%%bam}.BC.cpm.bs1.bedgraph -of bedgraph --binSize 1 --normalizeUsing CPM  --extendReads "
sbatch --time 12:00:00 --mem 4G   --wrap "bamCoverage --bam $file -o ${file%%bam}.BC.raw.bs1.bedgraph -of bedgraph --binSize 1  --extendReads "

done


mkdir seacr_no_ctrl_01

for file in *.bedgraph
do
bdg_1000=${file%%.bedgraph}
cmd_s_1000="SEACR_1.3.sh $file 0.01 non stringent seacr_no_ctrl_01/${bdg_1000}.0.01fdr"
echo $cmd_s_1000
sbatch -p epyc --mem 1G --time 01:00:00 --wrap "$cmd_s_1000"
done

for file in *.bedgraph
do
bdg_1000=${file%%.bedgraph}
cmd_s_1000="SEACR_1.3.sh $file 0.01 non relaxed seacr_no_ctrl_01/${bdg_1000}.rel.0.01fdr"
echo $cmd_s_1000
sbatch -p epyc --mem 1G --time 01:00:00 --wrap "$cmd_s_1000"
done
