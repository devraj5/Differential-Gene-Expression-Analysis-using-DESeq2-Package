# Differential-Gene-Expression-Analysis-using-DESeq2-Package
Analysed a RNA-seq data from a mammary adipose tissue in mouse model to study the role of Obesity in the poor breast cancer prognosis

Breast cancer (BC) is the most common tumor in women and represents the second cause of cancer-caused death after lung cancer 
In 2020, over 2.3 million new BC cases were estimated worldwide.
It has been recognized that the tumor stroma has an important role in cancer progression and dissemination. 
It provides nutrients, energy, and growth factors to comply with the metabolic requirements of cancer cells.
The predominant cell population in the breast is the adipocyte and its role in BC progression has only recently received significant attention.
Such adipocytes are termed cancer-associated adipocytes (CAAs) which have altered biological properties compared with the normal state
Within the breast cancer microenvironment, BC cells induce the production of endocrine and paracrine signaling mediators such as adipokines and metabolic substrates like bioactive lipids by CAAs. 
These in turn drive increased growth and invasion of tumor cells along with the development of resistance to chemotherapy and radiotherapy. Although, clear reciprocal interactions between adipocytes and BC cells are established, yet most of these studies have not considered the obesity context.
Much epidemiological evidence has demonstrated  that survival is lower in obese women, independent of menopausal status
Hence I had a research question that is- whether the influence of adipocytes on BC behavior is altered in obese patients compared with lean ones and could these significant differences help explain the observed epidemiological evidence

DATASET:

The gene expression data used in this research is from Gene Expression Omnibus(GEO) with the accession number – GSE201316
The dataset contained RNA expression profile of adipose tissues from 22 Female mice. Mice were fed either a regular chow diet(CD) or a high-fat diet (HFD) for 16 weeks and were injected with E0771 cancer cells resuspended in PBS mixed with Matrigel into their left inguinal/mammary fat pad (tumor-associated inguinal adipose tissue) while their right inguinal fat pad (control inguinal adipose tissue) received PBS-Matrigel solution alone.Tumor-associated (T) and control inguinal white adipose tissue (C) were collected two weeks after cancer cell injections 
RNA was harvested and sequenced using Illumina NextSeq 500.


Deseq2 Object:

There were two condtions and two genotypes:
Conditions/Treatments: T and C
Genotypes: HFD and CD

Running the Deseq2 object: dds <- DESeqDataSetFromMatrix(countData = main_data, colData = meta_data, design = ~genotype + condition + genotype:condition)
gave the result:
"Intercept"              "genotype_HFD_vs_CD"     "condition_T_vs_C"       "genotypeHFD.conditionT"

To observe the effect of treatment (tumor co-culture) in the obsese samples: Obese Tumor vs Obese control 
the contrast choose was:
l2f <- results(dds, lsit (c("condition_T_vs_C","genotypeHFD.conditionT")), alpha = 0.05)


To observe the effect of treatment (tumor co-culture) in the lean samples: lean Tumor vs lean control 
the contrast choose was:
l2f <- results(dds, contrast = c("condition","T", "C"), alpha = 0.05)

To observe the difference between lean and obsese with treatment (tumor co-culture): Obese Tumor vs lean tumor 
the contrast choose was:
l2f <- results(dds, list(c(  "genotype_HFD_vs_CD","T",  "genotypeHFD.conditionT"), alpha = 0.05)

To observe the difference between lean and obsese without treatment (control co-culture): Obese Control vs lean control 
the contrast choose was:
l2f <- results(dds, contrast =c("genotype_HFD_vs_CD","T",  "HFD", "CD")), alpha = 0.05)

And the final, the one I wanted to see. I wanted to select only those genes that were differentially expressed among all genotypes (HFD_T vs HFD_C vs CD_T vs CD_C) or lets say the difference response of treatment across genotype.  
the contrast choose was:
l2f <- results(dds, name="genotypeHFD.conditionT", alpha = 0.05)
