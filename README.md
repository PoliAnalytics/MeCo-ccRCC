ASSESSING THE ROLE OF MECHANOTRANSDUCTION IN ccRCC

---------------------------------------------------------------------------------------------------
ABOUT US 🌍
---------------------------------------------------------------------------------------------------
Hi 👋, we are Carlo Manenti👨🏻‍🔬 (Doc in Biotchenology), Paola Maragno👩🏼‍🔬(Doc in Biotechnology), Gabriele Marchi👨🏻‍🔬 (Doc in Biology) and Alberto Pettenella👨🏻‍🔬 (Doc in Bioengineering).
Our study focuses on  assessing  the role of mechanotransduction in clear cell Renal Cell Carcinoma (ccRCC).


---------------------------------------------------------------------------------------------------
INTRODUCTION 📚
---------------------------------------------------------------------------------------------------
ccRCC is the primary histological subtype of renal cancer. 
In particular, in our study, we focused on the role of cancer associated fibroblasts (CAFs), because CAFs are strongly involved in modifying the extracellular matrix (ECM) in the ccRCC. Changes in the ECM can lead to changes of gene expression in the cancer cells through mechanotransduction.

Mechanotransduction is a phenomenon through which mechanical cues, as the topology and stiffness of the matrix, are sense by the cell and translated into adaptive responses, as a switch in gene expression. 

In literature we found a way to quantify the effect of mechanotranduction on gene expression, the Mechanical Conditioning score (MeCo). In our study we adapt the MeCo for ccRCC and we are going to use it to asse the role of mechanical cues in the tumoral environment. 



---------------------------------------------------------------------------------------------------
 AIM OF THE STUDY 💡
---------------------------------------------------------------------------------------------------
The aim of the study is thus to firstly assess the power of MeCo to evaluate the condition of ccRCC patients. 

In this project we want to study mechanotransduction as the means to achieve a more effective therapy for the patient: consequently this analysis might interest pharmaceutical companies that want to have more more insight in the biology of ccRCC to develop treatments that target the key molecules and processes involved in it. 

Eventually we aim to develop a methodology for the analysis of tumour samples that could be used by clinicians in a versatile way: not only to study this type of tumour but ideally to be applied to every kind of cancer.



---------------------------------------------------------------------------------------------------
DATA PRE-PROCESSING 👩🏼‍💻🧑🏻‍💻
---------------------------------------------------------------------------------------------------
We used Gene Expression Omnibus to find gene expression data of CAFs cultured on soft and stiff substrates. For each condition we have 3 biological replicates. With DESeq2 we performed Differential gene expression analysis to see which genes are interested by mechanotransduction. To select differential expressed genes we used a Log2FoldChange of |3| and a p.adjusted value < 0.001.For FDR we used Benjamini-Hochberg and we removed all the lines that did not have more than 9 reads for sample. 

We needed to refine even further the MeCo so we relay on pathway analysis with Metascape. We selected the 5 most impactful pathways in tumor progression and compute a different MeCo scores for each of those pathway. 

So we obtain five refined MeCo scores that provide meaningful insights on tumoral pathways. 


---------------------------------------------------------------------------------------------------
MAIN DATASET 💾
---------------------------------------------------------------------------------------------------
For our dataset we use the TGCA-KIRC project, which is a Renal cancer data set. 
We selected only the patients that has both phenotypic and gene expression data, ending up with over 500 patients. 

For each we have phenotypic data, like: age, sex, stage and survival data. 

We also have gene expression data for most of the genes and so we were able to compute the refined MeCo scores for each patient. In the end we calculate a ECM MeCo, Proliferation MeCo, Antitumoral mechanisms MeCo, Inflammation MeCo and Chemotaxis MeCo. 



---------------------------------------------------------------------------------------------------
EXAMPLE CASE: ECM MeCO refined 🔎
---------------------------------------------------------------------------------------------------
ECM MeCO give us an idea of the ECM stiffness and shaping. The higher the score, the closer the  gene expression patter is to a stiff substrate condition. 

We can see that in older patients (60 =< years) we have a higher median than young (years < 60), as expected given the ageing process in older people.

But with a higher stage a lower value is observed. This is quite strage. Given the typical progression of a tumor we would expect a stiffer substrate with a higher stage. But in ccRenal Cancer this is the opposite. So the ECM MeCo seems able to follow nicely also this trend.



---------------------------------------------------------------------------------------------------
 RESEARCH PLAN 🚀
---------------------------------------------------------------------------------------------------
In the next few weeks we will use the MeCo scores to asses the role of mechanotransudction via hypothesis testing, survival analysis and other powerful prediction models






---------------------------------------------------------------------------------------------------
LAST. FINAL NOTES 👋
---------------------------------------------------------------------------------------------------
We will soon add also notes of the articles that we used as biological base for this research
Thank you for your time 😊
