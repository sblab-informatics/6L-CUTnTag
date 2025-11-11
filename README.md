# 6L-CUTnTag

Data processing codes for paper 6L-CUT&Tag

Requirements:

	•	demuxFQ (3.1.0)
	•	FastQC (0.11.8)
	•	MultiQC (1.18)
	•	Picard (2.20.3)
	•	samtools (1.18) 
	•	bedtools (2.31.0)
	•	bwa
	•	SEACR (1.3)
	•	deepTools (3.5.5)

The self-attention model was constructed with PyTorch (2.5.1) on CUDA (12.4). Other Python packages used in the model training, evaluation and testing includes:
- NumPy (1.26.4)
- SciPy (1.13.1)
- scikit-learn (1.5.2)
- torchsnooper (0.8)
- pickle (4.0)
	
Prediction:
```shell
python predict.py --input [your data path]
```
