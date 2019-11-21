There are five major sections in https://github.com/oxwang/fda_scRNA-seq

# Section 1: Preprocessing

In preprocessing section, we provided script code to process bulk RNA-Seq and single cell RNA-Seq (scRNA-Seq) fastq data.
For bulk RNA-Seq, one pipeline was provided. For scRNA-Seq, three pipelines were provided for 10x, C1_LLU, C1_FDA_HT, and
iCELL8 data, respectively.


# Section 2: Pipeline_comp

In pipeline_comp section, we provided R code to generate Figure 2 in our manuscript. The code is to compare the differences
between any two pipelines for scRNA-Seq data.


# Section 3: Normalization

In normalization section, we provided R code to generate sihouette values used for Figure 3 in our manuscript. The code is to
generate sihouette values for the two gene count matrices with different sequencing depth 10k and 100k for seven datasets
(10x_LLU, 10x_NCI, 10x_NCI_M, C1_LLU, C1_FDA_HT, iCELL8_PE, and iCELL8_SE) of HCC1395 and HCC1395BL samples. For each datasets
in HCC1395 and HCC1395BL, there is an R code file.


# Section 4: Batch correction

In batch correction section, we provided R and python code to perform batch correction. The 'batch.R' represents the data loading
code for four scenarios and Tian's data, and the code of seven batch correction. 'scanorama' folder represesnts the python code
used for scanorama batch correction. 'bbknn.py' file represents the python code for bbknn batch correction. 'metrics_scores.R'
represents the code to calculate alignment score, kBET score, and silhouette value.


# Section 5: others

In 'others' section, we provided the code to perform benchmarking assessment in section 5 & 6 in manuscript.
