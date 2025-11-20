#!/bin/bash
#SBATCH --mem 70G
#SBATCH -t 03:00:00



bedtools intersect -abam /scratchc/sblab/dhir01/Project/Rafa/Revision/Paper_bams/6LCnTH3K4Me1rep2.genome.GRCm38_primary_assembly.dedup.bam  -b no_chr_primed.bed >tmp.bam

samtools view -h tmp.bam | grep -v 'chrM' | grep -v GL | grep -v 'JH' >6LCnTH3K4Me1rep2primed.sam
samtools view -b 6LCnTH3K4Me1rep2primed.sam >6LCnT_H3K4me1_primed.bam
samtools index 6LCnT_H3K4me1_primed.bam



bedtools intersect -abam SRR4453260.merged.bam  -b primed.bed >tmp.bam
samtools view -h tmp.bam | sed 's/chr//g' | grep -v GL | grep -v 'JH' >SRR4453260.primed.sam
samtools view -b SRR4453260.primed.sam >CruzMolina_H3K4me1_primed.bam
samtools index CruzMolina_H3K4me1_primed.bam

#########################################




bedtools intersect -abam /scratchc/sblab/dhir01/Project/Rafa/Revision/Paper_bams/6LCnTH3K4Me1rep2.genome.GRCm38_primary_assembly.dedup.bam  -b no_chr_active.bed >tmp.bam

samtools view -h tmp.bam | grep -v 'chrM' | grep -v GL | grep -v 'JH' >6LCnTH3K4Me1rep2actve.sam
samtools view -b 6LCnTH3K4Me1rep2actve.sam >6LCnT_H3K4me1_active.bam
samtools index 6LCnT_H3K4me1_active.bam


bedtools intersect -abam /scratchc/sblab/dhir01/Project/Rafa/Revision/Paper_bams/6LCnTH3K27acrep2.genome.GRCm38_primary_assembly.dedup.bam  -b no_chr_active.bed >tmp.bam

samtools view -h tmp.bam | grep -v 'chrM' | grep -v GL | grep -v 'JH' >6LCnTH3K27acrep2active.sam
samtools view -b 6LCnTH3K27acrep2active.sam >6LCnT_H3K27ac_active.bam
samtools index 6LCnT_H3K27ac_active.bam


bedtools intersect -abam SRR4453260.merged.bam  -b active.bed >tmp.bam
samtools view -h tmp.bam | sed 's/chr//g' | grep -v GL | grep -v 'JH' >SRR4453260.active.sam
samtools view -b SRR4453260.active.sam >CruzMolina_H3K4me1_active.bam
samtools index CruzMolina_H3K4me1_active.bam


bedtools intersect -abam SRR4453258merged.bam  -b active.bed >tmp.bam
samtools view -h tmp.bam | sed 's/chr//g' | grep -v GL | grep -v 'JH' >SRR4453258.active.sam
samtools view -b SRR4453258.active.sam >CruzMolina_H3K27ac_active.bam
samtools index CruzMolina_H3K27ac_active.bam

#########################################


bedtools intersect -abam /scratchc/sblab/dhir01/Project/Rafa/Revision/Paper_bams/6LCnTH3K4Me1rep2.genome.GRCm38_primary_assembly.dedup.bam  -b no_chr_poised.bed >tmp.bam

samtools view -h tmp.bam | grep -v 'chrM' | grep -v GL | grep -v 'JH' >6LCnTH3K4Me1rep2poised.sam
samtools view -b 6LCnTH3K4Me1rep2poised.sam >6LCnT_H3K4me1_poised.bam
samtools index 66LCnT_H3K4me1_poised.bam



bedtools intersect -abam /scratchc/sblab/dhir01/Project/Rafa/Revision/Paper_bams/6LCnTH3K27Me3rep2.genome.GRCm38_primary_assembly.dedup.bam  -b no_chr_poised.bed >tmp.bam

samtools view -h tmp.bam | grep -v 'chrM' | grep -v GL | grep -v 'JH' >6LCnTH3K27Me3rep2poised.sam
samtools view -b 6LCnTH3K27Me3rep2poised.sam >6LCnT_H3K27me3_poised.bam
samtools index 6LCnT_H3K27me3_poised.bam


bedtools intersect -abam SRR4453260.merged.bam  -b poised.bed >tmp.bam
samtools view -h tmp.bam | sed 's/chr//g' | grep -v GL | grep -v 'JH' >SRR4453260.poised.sam
samtools view -b SRR4453260.poised.sam >CruzMolina_H3K4me1_poised.bam
samtools index CruzMolina_H3K4me1_poised.bam


bedtools intersect -abam SRR4453259.merged.bam  -b poised.bed >tmp.bam
samtools view -h tmp.bam | sed 's/chr//g' | grep -v GL | grep -v 'JH' >SRR4453259.poised.sam
samtools view -b SRR4453259.poised.sam >CruzMolina_H3K27me3_poised.bam
samtools index CruzMolina_H3K27me3_poised.bam

#########################################

#multiBamSummary bins --bamfiles 6LCnT_H3K4me1_active.bam 6LCnT_H3K27ac_active.bam CruzMolina_H3K4me1_active.bam CruzMolina_H3K27ac_active.bam 6LCnT_H3K4me1_poised.bam 6LCnT_H3K27me3_poised.bam CruzMolina_H3K4me1_poised.bam CruzMolina_H3K27me3_poised.bam -o result1.npz


multiBamSummary bins --bamfiles 6LCnT_H3K4me1_primed.bam CruzMolina_H3K4me1_primed.bam 6LCnT_H3K4me1_active.bam 6LCnT_H3K27ac_active.bam CruzMolina_H3K4me1_active.bam CruzMolina_H3K27ac_active.bam 6LCnT_H3K4me1_poised.bam 6LCnT_H3K27me3_poised.bam CruzMolina_H3K4me1_poised.bam CruzMolina_H3K27me3_poised.bam -o all_results.npz

plotCorrelation -in all_results.npz -c pearson -p heatmap --skipZeros -o Cruz_molina_pearson_allplot.pdf --plotNumbers
echo "DONE"

echo "DONE"
multiBamSummary bins --bamfiles 6LCnT_H3K4me1_active.bam CruzMolina_H3K4me1_active.bam -o H3K4me1_active.npz --outRawCounts H3K4me1_active.tab
multiBamSummary bins --bamfiles 6LCnT_H3K27ac_active.bam CruzMolina_H3K27ac_active.bam -o H3K27ac_active.npz --outRawCounts H3K27ac_active.tab
multiBamSummary bins --bamfiles 6LCnT_H3K4me1_poised.bam CruzMolina_H3K4me1_poised.bam -o H3K4me1_poised.npz  --outRawCounts H3K4me1_poised.tab
multiBamSummary bins --bamfiles 6LCnT_H3K27me3_poised.bam CruzMolina_H3K27me3_poised.bam -o H3K4me1_poised.npz  --outRawCounts H3K4me1_poised.tab

