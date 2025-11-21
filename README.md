# 6B-CUTnTag

Data processing codes for paper 6B-CUT&Tag

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
	•	Bismark (v0.24.0)
	•	Macs2 (2.2.7.1)


The self-attention model was constructed with PyTorch (2.5.1) on CUDA (12.4). Other Python packages used in the model training, evaluation and testing includes:
- NumPy (1.26.4)
- SciPy (1.13.1)
- scikit-learn (1.5.2)
- torchsnooper (0.8)
- pickle (4.0)
	
Prediction:
```shell
cd selfAttn
python predict.py --input [your data path] --cpt [select a model from checkpoints]
```
- ```checkpoints/CnT6L_model_at_epoch22.cpt```: model trianed on the 6B-CUT&Tag data.
- ```checkpoints/WG6L_model_at_epoch20.cpt```: model trained on the whole-genome data.
