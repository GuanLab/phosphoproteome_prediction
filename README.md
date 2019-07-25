## Machine Learning Empowers Phosphoproteome Prediction in Cancers

This is the package of our winning algorithm in the 2017 NCI-CPTAC DREAM Proteogenomics Challenge. 

background: [Proteogenomics Challenge](https://www.synapse.org/#!Synapse:syn8228304/wiki/)

see also: [Hongyang Li and Yuanfang Guan's 1st Place Solution](https://www.synapse.org/#!Synapse:syn11522015/wiki/496744) 

Please contact (hyangl@umich.edu or gyuanfan@umich.edu) if you have any questions or suggestions.

![Figure1](figure/fig1.png?raw=true "Title")

---

## Installation
Git clone a copy of code:
```
git clone https://github.com/GuanLab/phosphoproteome_prediction.git
```
## Required dependencies

* [R](https://www.r-project.org/) (3.4.3)
* [python](https://www.python.org) (3.6.5)
* [numpy](http://www.numpy.org/) (1.13.3). It comes pre-packaged in Anaconda.
* [scikit-learn](http://scikit-learn.org) (0.19.0) A popular machine learning package. It can be installed by:
```
pip install -U scikit-learn
```
## Dataset

All the omic data are 2D matrices, where columns are cancer samples and rows are genes/proteins/phosphorylation sites. 
The proteomic and phosphoproteomic data originally came from [CPTAC-breast](https://cptac-data-portal.georgetown.edu/cptac/s/S015) and [CPTAC-ovary](https://cptac-data-portal.georgetown.edu/cptac/s/S020).
The genomic data originally came from [TCGA-breast](https://portal.gdc.cancer.gov/legacy-archive/search/f?filters=%7B%22op%22:%22and%22,%22content%22:%5B%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22cases.project.program.name%22,%22value%22:%5B%22TCGA%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22cases.project.project_id%22,%22value%22:%5B%22TCGA-BRCA%22%5D%7D%7D%5D%7D) and [TCGA-ovary](https://portal.gdc.cancer.gov/legacy-archive/search/f?filters=%7B%22op%22:%22and%22,%22content%22:%5B%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22cases.project.program.name%22,%22value%22:%5B%22TCGA%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22cases.project.project_id%22,%22value%22:%5B%22TCGA-OV%22%5D%7D%7D%5D%7D). 

During the challenge, we directly downloaded these data from the challenge [website](https://www.synapse.org/#!Synapse:syn8228304/files/). Unfortunatelly, this download link s no longer available for unregistered users. We therefore provided examples of dummy data in the directory data/raw/. 
* retrospective_breast_CNA_sort_common_gene_16884.txt
* retrospective_breast_phospho_sort_common_gene_31981.txt
* retrospective_breast_proteome_sort_common_gene_10005.txt
* retrospective_breast_RNA_sort_common_gene_15107.txt
* retrospective_ova_CNA_sort_common_gene_11859.txt
* retrospective_ova_JHU_proteome_sort_common_gene_7061.txt
* retrospective_ova_phospho_filtered.txt
* retrospective_ova_phospho_sort_common_gene_10057.txt
* retrospective_ova_PNNL_proteome_sort_common_gene_7061.txt 
* retrospective_ova_rna_seq_sort_common_gene_15121.txt

Then preprocess the data and generate 5-fold cross validatation using code in 
* data/trimmed_set
* data/normalization
* data/cv_set

## Four models

We have two sets of code in parallel
* prediction/breast
* prediction/ova

### 1. "proteome" model
This model directly approximates the phosphorylation level based on the corresponding parent protein level. 
```
prediction/breast/proteome/
```

### 2. "site-specific" model
This model considers the protein-protein interactions in regulating phosphorylation, in which all protein levels were used as features to make predictions. The base learner is random forest with maximum depth of 3 and 100 trees.
```
prediction/breast/individual/
```

### 3. "cross-tissue" model
Similar to the "site-specific" model, this model uses combined samples from breast and ovarian cancer samples.
```
prediction/breast/individual_transplant/
```

### 4. "multi-site" model
This model considers the associations between phosphorylation sites of the same parent protein.
```
prediction/breast/multisite/
```

### 5. ensemble model
Our final results are the ensemble of the 1-4 models mentioned above.
```
prediction/breast/final/
```

### 6. result analysis and figure preparation
```
analysis_sub3/
```
