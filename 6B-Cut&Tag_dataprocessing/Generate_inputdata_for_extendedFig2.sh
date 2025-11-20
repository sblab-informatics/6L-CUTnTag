#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 25G
#SBATCH -t 02:00:00

#bedtools map -a ../../mm10_1kb_bins.bed -b H3K4Me1.merged.signal.bg -c 4 -o mean| LC_COLLATE=C sort  | awk '$4>=0' >H3K4Me1.merged.signal_1kb_bins.bg
#bedtools map -a ../../mm10_1kb_bins.bed -b H3K4Me3.merged.signal.bg -c 4 -o mean| LC_COLLATE=C sort  | awk '$4>=0' >H3K4Me3.merged.signal_1kb_bins.bg
#bedtools map -a ../../mm10_1kb_bins.bed -b H3K27Me3.merged.signal.bg -c 4 -o mean| LC_COLLATE=C sort  | awk '$4>=0' >H3K27Me3.merged.signal_1kb_bins.bg
#bedtools map -a ../../mm10_1kb_bins.bed -b H3K27ac.merged.signal.bg -c 4 -o mean| LC_COLLATE=C sort  | awk '$4>=0' >H3K27ac.merged.signal_1kb_bins.bg

bedtools intersect -a H3K4Me1.merged.signal_1kb_bins.bg -b H3K4Me1B2.mc.20x.bg -wa -wb | awk '{OFS="\t"; print $5"_"$6"_"$7,$8,$4}' >tmp1
wc -l tmp1
cat header tmp1 >H3K4Me1_enrichment.20x.mc.txt




bedtools intersect -a H3K4Me3.merged.signal_1kb_bins.bg -b H3K4Me3B2.mc.20x.bg -wa -wb | awk '{OFS="\t"; print $5"_"$6"_"$7,$8,$4}' >tmp1
wc -l tmp1
cat header tmp1 >H3K4Me3_enrichment.20x.mc.txt



bedtools intersect -a H3K27Me3.merged.signal_1kb_bins.bg -b H3K27Me3B2.mc.20x.bg -wa -wb | awk '{OFS="\t"; print $5"_"$6"_"$7,$8,$4}' >tmp1
wc -l tmp1
cat header tmp1 >H3K27Me3_enrichment.20x.mc.txt




bedtools intersect -a H3K27ac.merged.signal_1kb_bins.bg -b H3K27acB2.mc.20x.bg -wa -wb | awk '{OFS="\t"; print $5"_"$6"_"$7,$8,$4}' >tmp1
wc -l tmp1
cat header tmp1 >H3K27ac_enrichment.20x.mc.txt



bedtools intersect -a H3K4Me1.merged.signal_1kb_bins.bg -b WG.mc.20x.bg -wa -wb | awk '{OFS="\t"; print $5"_"$6"_"$7,$8,$4}' >tmp1
wc -l tmp1
cat header tmp1 >WG_003.H3K4Me1_enrichment.20x.mc.txt




bedtools intersect -a H3K4Me3.merged.signal_1kb_bins.bg -b WG.mc.20x.bg -wa -wb | awk '{OFS="\t"; print $5"_"$6"_"$7,$8,$4}' >tmp1
wc -l tmp1
cat header tmp1 >WG_003.H3K4Me3_enrichment.20x.mc.txt



bedtools intersect -a H3K27Me3.merged.signal_1kb_bins.bg -b WG.mc.20x.bg -wa -wb | awk '{OFS="\t"; print $5"_"$6"_"$7,$8,$4}' >tmp1
wc -l tmp1
cat header tmp1 >WG_003.H3K27Me3_enrichment.20x.mc.txt


bedtools intersect -a H3K27ac.merged.signal_1kb_bins.bg -b WG.mc.20x.bg -wa -wb | awk '{OFS="\t"; print $5"_"$6"_"$7,$8,$4}' >tmp1
wc -l tmp1
cat header tmp1 >WG_003.H3K27ac_enrichment.20x.mc.txt


#########################################

bedtools intersect -a H3K4Me1.merged.signal_1kb_bins.bg -b H3K4Me1B2.hmc.20x.bg -wa -wb | awk '{OFS="\t"; print $5"_"$6"_"$7,$8,$4}' >tmp2
wc -l tmp2
cat header tmp2 >H3K4Me1_enrichment.20x.hmc.txt




bedtools intersect -a H3K4Me3.merged.signal_1kb_bins.bg -b H3K4Me3B2.hmc.20x.bg -wa -wb | awk '{OFS="\t"; print $5"_"$6"_"$7,$8,$4}' >tmp2
wc -l tmp2
cat header tmp2 >H3K4Me3_enrichment.20x.hmc.txt



bedtools intersect -a H3K27Me3.merged.signal_1kb_bins.bg -b H3K27Me3B2.hmc.20x.bg -wa -wb | awk '{OFS="\t"; print $5"_"$6"_"$7,$8,$4}' >tmp2
wc -l tmp2
cat header tmp2 >H3K27Me3_enrichment.20x.hmc.txt




bedtools intersect -a H3K27ac.merged.signal_1kb_bins.bg -b H3K27acB2.hmc.20x.bg -wa -wb | awk '{OFS="\t"; print $5"_"$6"_"$7,$8,$4}' >tmp2
wc -l tmp2
cat header tmp2 >H3K27ac_enrichment.20x.hmc.txt



bedtools intersect -a H3K4Me1.merged.signal_1kb_bins.bg -b WG.hmc.20x.bg -wa -wb | awk '{OFS="\t"; print $5"_"$6"_"$7,$8,$4}' >tmp2
wc -l tmp2
cat header tmp2 >WG_003.H3K4Me1_enrichment.20x.hmc.txt




bedtools intersect -a H3K4Me3.merged.signal_1kb_bins.bg -b WG.hmc.20x.bg -wa -wb | awk '{OFS="\t"; print $5"_"$6"_"$7,$8,$4}' >tmp2
wc -l tmp2
cat header tmp2 >WG_003.H3K4Me3_enrichment.20x.hmc.txt



bedtools intersect -a H3K27Me3.merged.signal_1kb_bins.bg -b WG.hmc.20x.bg -wa -wb | awk '{OFS="\t"; print $5"_"$6"_"$7,$8,$4}' >tmp2
wc -l tmp2
cat header tmp2 >WG_003.H3K27Me3_enrichment.20x.hmc.txt


bedtools intersect -a H3K27ac.merged.signal_1kb_bins.bg -b WG.hmc.20x.bg -wa -wb | awk '{OFS="\t"; print $5"_"$6"_"$7,$8,$4}' >tmp2
wc -l tmp2
cat header tmp2 >WG_003.H3K27ac_enrichment.20x.hmc.txt
