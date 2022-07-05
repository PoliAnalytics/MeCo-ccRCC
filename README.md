# ASSESSING THE ROLE OF MECHANOTRANSDUCTION IN ccRCC

---------------------------------------------------------------------------------------------------
ABOUT US 🌍
---------------------------------------------------------------------------------------------------
Hi 👋, we are Carlo Manenti👨🏻‍🔬 (Doc in Biotechnology), Paola Maragno👩🏼‍🔬(Doc in Biomolecular Science and Technology) and Alberto Pettenella👨🏻‍🔬 (Doc in Bioengineering).
Our study focuses on  assessing  the role of mechanotransduction in clear cell Renal Cell Carcinoma (ccRCC) via computing __mechanical conditioning (MeCo) scores__.


---------------------------------------------------------------------------------------------------
INTRODUCTION 📚
--------------------------------------------------------------------------------------------------- 
ccRCC is the primary histological subtype of renal cancer. 
In particular, in our study, we focused on the role of cancer associated fibroblasts (CAFs), because CAFs are involved in modifying the extracellular matrix (ECM) in ccRCC. Changes in the ECM can lead to changes of gene expression in the cancer cells through mechanotransduction.

Mechanotransduction is a phenomenon through which mechanical cues, as the topology and stiffness of the matrix, are sense by the cell and translated into adaptive responses, as a switch in gene expression. 

From the literature we found a way to quantify the effect of mechanotranduction on gene expression, the Mechanical Conditioning score (MeCo). 
Here we adapt the MeCo for ccRCC and we use it to assess the role of mechanical cues in the tumoral environment. 



---------------------------------------------------------------------------------------------------
 AIM OF THE STUDY 💡
---------------------------------------------------------------------------------------------------
The aim of this study is thus to firstly assess the power of MeCo to evaluate the condition of ccRCC patients. 

In this project we want to study mechanotransduction as the means to achieve a more effective therapy for the patient: consequently we hope that this analysis will be of interest for pharmaceutical companies wanting to gain an insight in the biology of ccRCC to develop treatments targeting the key molecules and processes involved in it. 

Eventually, provided that MeCo assessment yields consistent results, we wish to develop a methodology for the analysis of tumor samples that could be benefit clinicians by not being restricted only to ccRCC.


---------------------------------------------------------------------------------------------------
DATA PRE-PROCESSING 👩🏼‍💻🧑🏻‍💻🧑🏻‍💻
---------------------------------------------------------------------------------------------------
We used Gene Expression Omnibus to find gene expression data of CAFs cultured on soft and stiff substrates. For each condition we have 3 biological replicates. First of all we converted all the Ensembl ID to Gene ID (GRCh38) and discarded rRNA, mitochondrial and non reliable gene annotations. Than using edgeR we performed Differential Gene Expression Analysis to identify which genes are interested by mechanotransduction. To select differential expressed genes we used a Log2FoldChange of |1| and a p.adjusted value < 0.001. While to adjust for FDR, we used Benjamini-Hochberg. 

To refine even further the MeCo so we relay on pathway analysis with Metascape. 
We selected the 3 most impactful pathways, related to the selected genes, in tumor progression and computed MeCo scores for each of those pathways: MeCo Development, MeCo regulation (of biological process) and MeCo response (to stimulus). 

---------------------------------------------------------------------------------------------------
MAIN DATASET 💾
---------------------------------------------------------------------------------------------------
As the chosen dataset, we used the TGCA-KIRC project, which is a Renal cancer data set. 
We selected only the patients that have both phenotypic and gene expression data, ending up with over 500 patients. 

For each individual we kept as variables of interest related to phenotypic data: age, sex, tumor stage and survival data. 

MeCo scores were computed for each patient given the overlap of the previously identified genes and the transcriptomics data provided in the dataset.


---------------------------------------------------------------------------------------------------
EXPLORATIVE ANALYSIS: 'There can only be one' ⚡️
---------------------------------------------------------------------------------------------------
Using the Mann whitney U test we assed if the MeCo and its iterations (MeCo refined) were able to follow both physiological and pathological trends, like a higher stiffness due to aging and lower stiffness in higher tumoral stages due to a peculiarity of ccRCC. Eventually only the MeCo (general) and MeCo regulation (refined) were able to follow nicely those trend. 
So we keep them for the down stream analysis. 

---------------------------------------------------------------------------------------------------
 RESEARCH PLAN 🚀
---------------------------------------------------------------------------------------------------
- Survival Analysis ⚕️
- Multi class Logistic Regression 🏷
- Mixed Effect Model 🎛


---------------------------------------------------------------------------------------------------
- Survival Analysis ⚕️
---------------------------------------------------------------------------------------------------
In the survival analysis, both MeCo and MeCo refined were highlighted as protective factors. 

Cox PH model:
We defined a Cox PH model accounting for Age, Stage, and MeCo or MeCo refined. After stratifying 
for the stage the assumptions for the Cox PH model were satisfied and were obtain a reduction of 
around 2.4% in the risk of death with an increase of 0.01 of the MeCo| MeCo refined score. 

note: We evaluated the model using a 0.01 increase since both MeCo and MeCo refined scores range from -0.25, up to 0.50. 
So a unitary increase would mean a leap from one extreme to the other of the MeCo scores. 

Also, we used an exponential accelerated failure time model to assess the role of MeCo, 
verifying once again the results obtained by the Cox PH models. 

---------------------------------------------------------------------------------------------------
- Multi class Logistic Regression 🏷
---------------------------------------------------------------------------------------------------
We also used a multi-class logistic regression model to predict the tumor stage from the age, the status, and MeCo | MeCo refined. For this task, we obtain interesting results since the MeCo and MeCo refined were meaningful variables for stages 2, 3, and 4 with consideration of the baseline (stage 1). 


---------------------------------------------------------------------------------------------------
- Mixed Effect Model 🎛
---------------------------------------------------------------------------------------------------
Also, we used logistic regression to predict the status of the patients given the age, stage, and Meco | MeCo refined. To further assess the role played by the stage variable we applied a Mixed Effect Model with a random intercept to assess if and how a given stage weight on the final result. 

As expected different stages account for almost 30% of the variability of the final result. 

---------------------------------------------------------------------------------------------------
CONCLUSIONS🕵🏼‍♀️🕵🏻🕵🏻
---------------------------------------------------------------------------------------------------
In the end, mechanotransduction plays are relevant role in ccRCC accounting for both survival time and tumoral stage. 

Validation with 'wet' lab techniques is required but we are confident that this kind of analysis can be easily extended to other types of tumors and may give a new perspective on the problem. 

---------------------------------------------------------------------------------------------------
REFERENCES 📚
---------------------------------------------------------------------------------------------------
Bond KH, Chiba T, Wynne KPH, et al. , The Extracellular Matrix Environment of Clear Cell Renal Cell Carcinoma Determines Cancer Associated Fibroblast Growth. Cancers (Basel), 2021 

Watson, Adam W. et al., Breast tumor stiffness instructs bone metastasis via maintenance of mechanical conditioning, Cell Reports, Volume 35, Issue 13, 109293, 2021

Paradiso, F., Quintela, M., Lenna, S., Serpelloni, S., James, D., Caserta, S., Conlan, S., Francis, L. and Taraballi, F., Studying Activated Fibroblast Phenotypes and Fibrosis-Linked Mechanosensing Using 3D Biomimetic Models. Macromol. Biosci., 22: 2100450. 2022

