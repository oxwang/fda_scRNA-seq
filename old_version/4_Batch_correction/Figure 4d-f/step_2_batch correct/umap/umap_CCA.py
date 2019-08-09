import os
import scanpy.api as sc


#2dataset
test = sc.read_csv('./CCA_sean-pcs/cca_aligned.csv')
print(test)
sc.pp.neighbors(test)
sc.tl.umap(test)
f = open('./CCA_sean-pcs/umap.csv', 'w')
for i in test.obsm:
    f.write(','.join([str(j) for j in i[0]]) + '\n')
f.close()


#2dataset-tumor
test = sc.read_csv('./CCA_sean-pcs/cca_aligned.csv')
print(test)
sc.pp.neighbors(test)
sc.tl.umap(test)
f = open('./CCA_sean-pcs/umap.csv', 'w')
for i in test.obsm:
    f.write(','.join([str(j) for j in i[0]]) + '\n')
f.close()


#4dataset
test = sc.read_csv('./CCA_sean-pcs/cca_aligned.csv')
print(test)
sc.pp.neighbors(test)
sc.tl.umap(test)
f = open('./CCA_sean-pcs/umap.csv', 'w')
for i in test.obsm:
    f.write(','.join([str(j) for j in i[0]]) + '\n')
f.close()


#LLU_NCI
test = sc.read_csv('./CCA_sean-pcs/cca_aligned.csv')
print(test)
sc.pp.neighbors(test)
sc.tl.umap(test)
f = open('./CCA_sean-pcs/umap.csv', 'w')
for i in test.obsm:
    f.write(','.join([str(j) for j in i[0]]) + '\n')
f.close()
