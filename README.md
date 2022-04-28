# ASSESSING THE ROLE OF MECHANOTRANSDUCTION IN ccRCC

---------------------------------------------------------------------------------------------------
ABOUT US 🌍
---------------------------------------------------------------------------------------------------
Hi 👋, we are Carlo Manenti👨🏻‍🔬 (Doc in Biotchenology), Paola Maragno👩🏼‍🔬(Doc in Biomolecular Science and Technology), Gabriele Marchi👨🏻‍🔬 (Doc in Biology) and Alberto Pettenella👨🏻‍🔬 (Doc in Bioengineering).
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
The aim of the study is thus to firstly assess the power of MeCo to evaluate the condition of ccRCC patients. 

In this project we want to study mechanotransduction as the means to achieve a more effective therapy for the patient: consequently we hope that this analysis will be of interest for pharmaceutical companies wanting to gain an insight in the biology of ccRCC to develop treatments targeting the key molecules and processes involved in it. 

Eventually, provided that MeCo assessment yields consistent results, we wish to develop a methodology for the analysis of tumor samples that could be benefit clinicians by not being restricted to ccRCC only.


---------------------------------------------------------------------------------------------------
DATA PRE-PROCESSING 👩🏼‍💻🧑🏻‍💻
---------------------------------------------------------------------------------------------------
We used Gene Expression Omnibus to find gene expression data of CAFs cultured on soft and stiff substrates. For each condition we have 3 biological replicates. With DESeq2 we performed Differential Gene Expression Analysis to identify which genes are interested by mechanotransduction. To select differential expressed genes we used a Log2FoldChange of |3| and a p.adjusted value < 0.001. 
To control for FDR, we used Benjamini-Hochberg and we removed all the genes that did not have more than 9 reads for sample. 

To refine even further the MeCo so we relay on pathway analysis with Metascape. 
We selected the 5 most impactful pathways, related to the selected genes, in tumor progression and computed MeCo scores for each of those pathways: ECM MeCo, Proliferation MeCo, Antitumoral mechanisms MeCo, Inflammation MeCo and Chemotaxis MeCo. 

---------------------------------------------------------------------------------------------------
MAIN DATASET 💾
---------------------------------------------------------------------------------------------------
As the chosen dataset, we used the TGCA-KIRC project, which is a Renal cancer data set. 
We selected only the patients that have both phenotypic and gene expression data, ending up with over 500 patients. 

For each individual we kept as variables of interest related to phenotypic data: age, sex, tumor stage and survival data. 

MeCo scores were computed for each patient given the overlap of the previously identified genes and the transcriptomics data provided in the dataset.


---------------------------------------------------------------------------------------------------
EXAMPLE CASE: ECM MeCO refined 🔎
---------------------------------------------------------------------------------------------------
ECM MeCo give us an idea of the ECM stiffness and shaping. The higher the score, the closer the gene expression pattern is to a stiff substrate condition. 

We can see that in older patients (60 =< years) we have a higher median than in young ones (years < 60), as expected given the aging process.

But with a higher stage a lower value is observed. Given the typical progression of a tumor we would expect a stiffer substrate with a higher stage. But in ccRenal Cancer this is the opposite. So the ECM MeCo seems able to follow nicely also this trend.


---------------------------------------------------------------------------------------------------
 RESEARCH PLAN 🚀
---------------------------------------------------------------------------------------------------
In the next few weeks we will use the MeCo scores to asses the role of mechanotransduction via hypothesis testing, survival analysis and other powerful prediction models



---------------------------------------------------------------------------------------------------
REFERENCES 📚
---------------------------------------------------------------------------------------------------
Bond KH, Chiba T, Wynne KPH, et al. , The Extracellular Matrix Environment of Clear Cell Renal Cell Carcinoma Determines Cancer Associated Fibroblast Growth. Cancers (Basel), 2021 

Watson, Adam W. et al., Breast tumor stiffness instructs bone metastasis via maintenance of mechanical conditioning, Cell Reports, Volume 35, Issue 13, 109293, 2021

Paradiso, F., Quintela, M., Lenna, S., Serpelloni, S., James, D., Caserta, S., Conlan, S., Francis, L. and Taraballi, F., Studying Activated Fibroblast Phenotypes and Fibrosis-Linked Mechanosensing Using 3D Biomimetic Models. Macromol. Biosci., 22: 2100450. 2022

