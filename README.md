# NEP100 : Methodological Distinction from Prior Work

This repository supports the study:  
**"Multi-omics analysis constructs a novel neuroendocrine prostate cancer classifier and classification system"** (Sci Rep 2025).  

## Key Innovations vs. [Theranostics 2024, 14(3):1065-1080]
While leveraging **field-standard pipelines** (e.g., Seurat, DESeq2, random forest), our work introduces:

### 1. Methodological Advancements:(It conforms to the conventional pattern of methodological progress)
- **Unprecedented scale**:  
  Integrated **14 scRNA-seq datasets** + **19 bulk cohorts (n=3,000+)** → *Largest NEPC multi-omics pipeline to date*
- **Novel analytical frameworks**:  
  - A more reasonable feature selection method  (5 methods)
  - The quantitative model of tumor heterogeneity was proposed for the first time (The model based on Shannon entropy : THEnt)
  - A more reliable model establishment based on AUC comparison  (Weighted sum based on the importance ranking of RF features)

### 2. New Discovery of NEPC Classification system Based on NEP100: (It conforms to the cumulative characteristics of academic development)
- **First 4-subtype NEPC classification** (`VR_O`, `Prol_N`, `Prol_P`, `EMT_Y`) with:
  - Comprehensive gene set enrichment and variation analysis
  - Master regulator identification (CaCTS)
  - Pseudotime trajectory mapping of subtype evolution  
  - Drug sensitivity profiling (oncoPredict) → **VR_O chemotherapy resistance**  
- **AMIGO2 as a novel target**:  
  - Validated in **7 independent cohorts + 10-patient IHC**  
  - Mechanistically linked to ADT resistance  

### 3. Independent Validation Tier
- **Public cohorts** The reproducibility of key metrics (e.g., The key genes are highly enriched in the tumor region of NPEC) using identical public cohorts (e.g., GSE230282 in Figure 4) is expected and desirable in computational biology, as it confirms the robustness of analytical frameworks.
- **Added clinical translation layer** absent in prior computational studies.
