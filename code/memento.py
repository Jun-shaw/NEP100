import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import matplotlib.pyplot as plt
import pandas as pd
from scipy.sparse import csr_matrix

import scanpy as sc
import memento

data_path = './'
adata = sc.read(data_path + 'seurat.harmony.h5ad')
print(adata)
adata.X = csr_matrix(adata.X)

adata.obs['subtype'] = adata.obs['celltype12'].apply(lambda x: 1 if x == 1 else 0)

print(adata.obs[['celltype12','subtype']].sample(5))

result_1d = memento.binary_test_1d(
   adata=adata,
    capture_rate=0.25,
    treatment_col='subtype',
    num_cpus=12,
    num_boot=2000)

print(result_1d.query('de_coef > 0').sort_values('de_pval').head(10))

df = pd.DataFrame(result_1d)
df.to_csv('output_test.csv', index=False)