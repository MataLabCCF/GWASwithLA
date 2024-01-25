import os

imputationFolder = "/home/peixott/beegfs/Analysis/LocalAncestry/ImputationTest/RephasedImputation"
outFolder = "/home/peixott/beegfs/Analysis/LocalAncestry/TRACTOR"
plink2 = "/home/peixott/beegfs/Programs/plink2"

fileToMerge = open(f"{outFolder}/AllPGEN.txt", "w")

for chrom in range(1,23):
    os.system(f"{plink2} --vcf {imputationFolder}/chr{chrom}.dose.vcf.gz --make-pgen --out {outFolder}/typedChr{chrom} --require-no-info \"IMPUTED\"")
    fileToMerge.write(f"{outFolder}/typedChr{chrom}\n")

fileToMerge.close()
os.system(f"{plink2} --pmerge-list {outFolder}/AllPGEN.txt --make-pgen --out {outFolder}/typedAlLChrom")
os.system(f"{plink2} --pfile {outFolder}/typedAlLChrom --indep-pairwise 200 50 0.2 --out {outFolder}/toPruneLD")
os.system(f"{plink2} --pfile {outFolder}/typedAlLChrom --keep toPruneLD.prune.in --make-pgen --out {outFolder}/typedAllNoLD")
os.system(f"{plink2} --pfile {outFolder}/typedAlNoLD --pca 10 --out {outFolder}/PCA")


