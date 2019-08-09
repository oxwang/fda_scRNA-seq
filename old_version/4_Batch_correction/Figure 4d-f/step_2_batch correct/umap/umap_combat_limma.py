import os
import scanpy.api as sc


#2dataset
test = sc.read_csv('./ComBat/combat.csv')
test = test.transpose()
print(test)
sc.pp.pca(test, n_comps=11)
sc.pp.neighbors(test)
sc.tl.umap(test)
print(test.obsm['X_umap'])
f = open('./ComBat/umap-11pcs.csv', 'w')
for i in test.obsm['X_umap']:
    f.write(','.join([str(j) for j in i]) + '\n')
f.close()

test = sc.read_csv('./limma/limma.csv')
test = test.transpose()
print(test)
sc.pp.pca(test, n_comps=10)
sc.pp.neighbors(test)
sc.tl.umap(test)
print(test.obsm['X_umap'])
f = open('./limma/umap-10pcs.csv', 'w')
for i in test.obsm['X_umap']:
    f.write(','.join([str(j) for j in i]) + '\n')
f.close()


#2dataset-tumor
test = sc.read_csv('./ComBat/combat.csv')
test = test.transpose()
print(test)
sc.pp.pca(test, n_comps=8)
sc.pp.neighbors(test)
sc.tl.umap(test)
print(test.obsm['X_umap'])
f = open('./ComBat/umap-8pcs.csv', 'w')
for i in test.obsm['X_umap']:
    f.write(','.join([str(j) for j in i]) + '\n')
f.close()

test = sc.read_csv('./limma/limma.csv')
test = test.transpose()
print(test)
sc.pp.pca(test, n_comps=7)
sc.pp.neighbors(test)
sc.tl.umap(test)
print(test.obsm['X_umap'])
f = open('./limma/umap-7pcs.csv', 'w')
for i in test.obsm['X_umap']:
    f.write(','.join([str(j) for j in i]) + '\n')
f.close()


#4dataset
test = sc.read_csv('./ComBat/combat.csv')
test = test.transpose()
print(test)
sc.pp.pca(test, n_comps=10)
sc.pp.neighbors(test)
sc.tl.umap(test)
print(test.obsm['X_umap'])
f = open('./ComBat/umap-10pcs.csv', 'w')
for i in test.obsm['X_umap']:
    f.write(','.join([str(j) for j in i]) + '\n')
f.close()

test = sc.read_csv('./limma/limma.csv')
test = test.transpose()
print(test)
sc.pp.pca(test, n_comps=7)
sc.pp.neighbors(test)
sc.tl.umap(test)
print(test.obsm['X_umap'])
f = open('./limma/umap-7pcs.csv', 'w')
for i in test.obsm['X_umap']:
    f.write(','.join([str(j) for j in i]) + '\n')
f.close()


#LLU_NCI
test = sc.read_csv('./ComBat/combat.csv')
test = test.transpose()
print(test)
sc.pp.pca(test, n_comps=11)
sc.pp.neighbors(test)
sc.tl.umap(test)
print(test.obsm['X_umap'])
f = open('./ComBat/umap-11pcs.csv', 'w')
for i in test.obsm['X_umap']:
    f.write(','.join([str(j) for j in i]) + '\n')
f.close()

test = sc.read_csv('./limma/limma.csv')
test = test.transpose()
print(test)
sc.pp.pca(test, n_comps=7)
sc.pp.neighbors(test)
sc.tl.umap(test)
print(test.obsm['X_umap'])
f = open('./limma/umap-7pcs.csv', 'w')
for i in test.obsm['X_umap']:
    f.write(','.join([str(j) for j in i]) + '\n')
f.close()
