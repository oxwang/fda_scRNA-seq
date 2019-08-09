import os
import scanpy.api as sc


#2dataset
test = sc.read_csv('./mnnCorrect/logCounts.csv')
test = test.transpose()
print(test)
sc.pp.pca(test, n_comps=18)
sc.pp.neighbors(test)
sc.tl.umap(test)
print(test.obsm['X_umap'])
f = open('./mnnCorrect/uncorrected_umap-18pcs.csv', 'w')
for i in test.obsm['X_umap']:
    f.write(','.join([str(j) for j in i]) + '\n')
f.close()


test = sc.read_csv('./mnnCorrect/mnnCounts.csv')
test = test.transpose()
print(test)
sc.pp.pca(test, n_comps=19)
sc.pp.neighbors(test)
sc.tl.umap(test)
print(test.obsm['X_umap'])
f = open('./mnnCorrect/umap-19pcs.csv', 'w')
for i in test.obsm['X_umap']:
    f.write(','.join([str(j) for j in i]) + '\n')
f.close()


#2dataset-tumor
test = sc.read_csv('./mnnCorrect/logCounts.csv')
test = test.transpose()
print(test)
sc.pp.pca(test, n_comps=24)
sc.pp.neighbors(test)
sc.tl.umap(test)
print(test.obsm['X_umap'])
f = open('./mnnCorrect/uncorrected_umap-24pcs.csv', 'w')
for i in test.obsm['X_umap']:
    f.write(','.join([str(j) for j in i]) + '\n')
f.close()


test = sc.read_csv('./mnnCorrect/mnnCounts.csv')
test = test.transpose()
print(test)
sc.pp.pca(test, n_comps=26)
sc.pp.neighbors(test)
sc.tl.umap(test)
print(test.obsm['X_umap'])
f = open('./mnnCorrect/umap-26pcs.csv', 'w')
for i in test.obsm['X_umap']:
    f.write(','.join([str(j) for j in i]) + '\n')
f.close()


#4dataset
test = sc.read_csv('./mnnCorrect/logCounts.csv')
test = test.transpose()
print(test)
sc.pp.pca(test, n_comps=30)
sc.pp.neighbors(test)
sc.tl.umap(test)
print(test.obsm['X_umap'])
f = open('./mnnCorrect/uncorrected_umap-30pcs.csv', 'w')
for i in test.obsm['X_umap']:
    f.write(','.join([str(j) for j in i]) + '\n')
f.close()


test = sc.read_csv('./mnnCorrect/mnnCounts.csv')
test = test.transpose()
print(test)
sc.pp.pca(test, n_comps=34)
sc.pp.neighbors(test)
sc.tl.umap(test)
print(test.obsm['X_umap'])
f = open('./mnnCorrect/umap-34pcs.csv', 'w')
for i in test.obsm['X_umap']:
    f.write(','.join([str(j) for j in i]) + '\n')
f.close()


#LLU_NCI
test = sc.read_csv('./mnnCorrect/mnnCounts.csv')
test = test.transpose()
print(test)
sc.pp.pca(test, n_comps=39)
sc.pp.neighbors(test)
sc.tl.umap(test)
print(test.obsm['X_umap'])
f = open('./mnnCorrect/umap-39pcs.csv', 'w')
for i in test.obsm['X_umap']:
    f.write(','.join([str(j) for j in i]) + '\n')
f.close()

test = sc.read_csv('./mnnCorrect/logCounts.csv')
test = test.transpose()
print(test)
sc.pp.pca(test, n_comps=38)
sc.pp.neighbors(test)
sc.tl.umap(test)
print(test.obsm['X_umap'])
f = open('./mnnCorrect/uncorrected_umap-38pcs.csv', 'w')
for i in test.obsm['X_umap']:
    f.write(','.join([str(j) for j in i]) + '\n')
f.close()
