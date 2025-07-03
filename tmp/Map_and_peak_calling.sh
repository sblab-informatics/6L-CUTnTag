for file in *.r1.*.fastq
do
fq1=$file
fq2=${fq1/r1/r2}
ls $fq1
ls $fq2
echo "----"
sbatch --time 01:00:00 -o %j.out -e %j.err --mem 16000 --wrap "cutadapt -q 20  -o trimmed2/${fq1%%.fastq}.trimmed.fq.gz -p trimmed2/${fq2%%.fastq}.trimmed.fq.gz $fq1 $fq2 "
done

mkdir aligned
g='/Users/dhir01/reference_genomes/mm10.fa'
w='/Users/dhir01/reference_genomes/mm10.whitelist.sorted.bed'
  
#aln with bwa
for file in *.r1*trimmed.fq.gz 
  do
  f1=$file
  f2=${file/r1/r2}
  sbatch --time 12:00:00  -p epyc --mem 16G --wrap "bwa mem -M $g $f1 $f2  | samtools view -Sb -F780 -q 10 -L $w - | samtools sort -@ 8 -   > aligned/${f1%%.trimmed.fq.gz}.mm10.sort.bam"
done

for file in *.r1.00.mm10.sort.bam
do
sample=${file%%.r1.00.mm10.sort.bam}
sbatch --time 01:00:00 -p epyc --mem 8G --wrap "samtools merge $sample.merged.bam ${file%%.r1.00.mm10.sort.bam}.r1.*.mm10.sort.bam"
echo "-------"
done


for bam in *.bam
do
sbatch -p epyc --time 12:00:00 --mem 8G --wrap "java -Xmx7g -jar /home/simeon01/applications/picard-2.20.3.jar MarkDuplicates INPUT=$bam OUTPUT=${bam%%.bam}.markduplicates.bam REMOVE_DUPLICATES=true  AS=true METRICS_FILE=${bam%%.bam}.markduplicates_metrics.txt | samtools index ${bam%%.bam}.markduplicates.bam"
done


for file in *merged*.bam
do
#sbatch --time 12:00:00 --mem 4G  -p epyc --wrap "bamCoverage --bam $file -o ${file%%bam}.s150.cpm.bs10.bw --binSize 10 --normalizeUsing CPM  --extendReads --maxFragmentLength 150 --smoothLength 9 "
sbatch --time 02:00:00 --mem 4G  -p epyc --wrap "bamCoverage --bam $file -o ${file%%bam}.150.cpm.bs10.bw --binSize 10 --normalizeUsing CPM  --extendReads --maxFragmentLength 150"
sbatch --time 02:00:00 --mem 4G  -p epyc --wrap "bamCoverage --bam $file -o ${file%%bam}.120.cpm.bs10.bw --binSize 10 --normalizeUsing CPM  --extendReads --maxFragmentLength 120"
sbatch --time 02:00:00 --mem 4G  -p epyc --wrap "bamCoverage --bam $file -o ${file%%bam}.130.cpm.bs10.bw --binSize 10 --normalizeUsing CPM  --extendReads --maxFragmentLength 130"

#sbatch --time 12:00:00 --mem 4G  -p epyc --wrap "bamCoverage --bam $file -o ${file%%bam}.s350.cpm.bs10.bw --binSize 10 --normalizeUsing CPM  --extendReads --maxFragmentLength 350 --smoothLength 9 "
sbatch --time 02:00:00 --mem 4G -p epyc --wrap "bamCoverage --bam $file -o ${file%%bam}.350.cpm.bs10.bw --binSize 10 --normalizeUsing CPM  --extendReads --maxFragmentLength 350"

#sbatch --time 12:00:00 --mem 4G  -p epyc --wrap "bamCoverage --bam $file -o ${file%%bam}.s600.cpm.bs10.bw --binSize 10 --normalizeUsing CPM  --extendReads --maxFragmentLength 600 --smoothLength 9 "
#sbatch --time 02:00:00 --mem 4G  -p epyc --wrap "bamCoverage --bam $file -o ${file%%bam}.600.cpm.bs10.bw --binSize 10 --normalizeUsing CPM  --extendReads --maxFragmentLength 600"

#sbatch --time 12:00:00 --mem 4G  -p epyc --wrap "bamCoverage --bam $file -o ${file%%bam}.s250.cpm.bs10.bw --binSize 10 --normalizeUsing CPM  --extendReads --maxFragmentLength 250 --smoothLength 9 "
sbatch --time 02:00:00 --mem 4G  -p epyc --wrap "bamCoverage --bam $file -o ${file%%bam}.250.cpm.bs10.bw --binSize 10 --normalizeUsing CPM  --extendReads --maxFragmentLength 250"
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
