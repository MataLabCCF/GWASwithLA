import argparse
import os

def createSets(eigenvec, setSize, outputFolder, outputName, logFile):

    #Code from https://stackoverflow.com/questions/5452576/k-means-algorithm-variation-with-equal-cluster-size?rq=1
    from sklearn.cluster import KMeans
    from scipy.spatial.distance import cdist
    from scipy.optimize import linear_sum_assignment
    import numpy as np
    import pandas as pd

    X = pd.read_table(eigenvec, header=None, names=['IID', 'FID', 'PC1', 'PC2'], sep=" ", index_col="IID",
                      usecols=['IID', 'PC1', 'PC2'])
    n_clusters = int(np.ceil(len(X)/setSize))
    kmeans = KMeans(n_clusters)
    kmeans.fit(X)
    centers = kmeans.cluster_centers_
    centers = centers.reshape(-1, 1, X.shape[-1]).repeat(setSize, 1).reshape(-1, X.shape[-1])
    distance_matrix = cdist(X, centers)
    clusters = linear_sum_assignment(distance_matrix)[1]//setSize

    X['ClusterID'] = clusters
    X.to_csv(f"{outputFolder}/{outputName}_PCA_Cluster.csv")

    return f"{outputFolder}/{outputName}_PCA_Cluster.csv"

def phaseWithEagle(allSamplesVCF, referencePhase, outputFolder, begin, end, threads, geneticMap, eagle, logFile, run = True):
    for i in range(begin, end + 1):
        samplesWithChr = allSamplesVCF.replace('*', str(i))
        vcfRefWithChr = referencePhase.replace('*', str(i))
        command = f'{eagle} --vcfTarget {samplesWithChr} --outPrefix {outputFolder}/Merged_Phased_{i} --numThreads {threads}' \
                  f' --geneticMapFile {geneticMap} --keepMissingPloidyX'
        if referencePhase != '':
            command = f'{command} --vcfRef {vcfRefWithChr} --allowRefAltSwap'
        execute(command, logFile, run)
    return f'{outputFolder}/Merged_Phased_{i}.vcf.gz'


def calculatePCAWithPlink(vcf, plink, folder, begin, end, outputName, logFile, run = True):
    logFile.write(f"\n\n============================= PLINK PCA =============================\n\n")

    fileMerge = open(f'{folder}/ListToConcat.txt', 'w')
    for i in range(begin, end + 1):
        execute(f"{plink} --vcf {vcf} --make-bed --out {folder}/{outputName}_chr{i}", logFile, run)
        fileMerge.write(f"{folder}/{outputName}_chr{i}\n")
    fileMerge.close()

    execute(f"{plink} --merge-list {folder}/ListToConcat.txt --make-bed --out {folder}/{outputName}_All", logFile, run)
    execute(f"{plink} --bfile {folder}/{outputName}_All --indep-pairwise 200 50 0.2 --out {folder}/{outputName}_LD",
            logFile, run)
    execute(f"{plink} --bfile {folder}/{outputName}_All --extract {folder}/{outputName}_LD.prune.in --make-bed --out "
            f"{folder}/{outputName}_LD_Prunned", logFile, run)
    execute(f"{plink} --bfile {folder}/{outputName}_LD_Prunned --pca 2 --out {folder}/{outputName}_PCA", logFile, run)

    return f'{folder}/{outputName}_PCA.eigenvec'


def mergeVCFs(vcf1, vcf2, bcftools, folder, begin, end, thread, logFile, run = True):
    logFile.write(f"\n\n============================= Isec and Merge =============================\n\n")
    for i in range(begin, end + 1):
        vcf1WithChr = vcf1.replace("*", str(i))
        vcf2WithChr = vcf2.replace("*", str(i))

        execute(f"{bcftools} isec -n=2 {vcf1WithChr} {vcf2WithChr} -o {folder}/CommonVariants_chr{i} --threads {thread}", logFile, run)
        execute(f"{bcftools} merge -R {folder}/CommonVariants_chr{i} {vcf1WithChr} {vcf2WithChr} -Oz "
                f"-o {folder}/Merged_chr{i}.vcf.gz --threads {thread}", logFile, run)
        execute(f"{bcftools} index {folder}/Merged_chr{i}.vcf.gz --threads {thread}", logFile, run)
    return f"{folder}/Merged_chr*.vcf.gz"


def normAndBiallelic(vcfFile, folder, bcftools, begin, end, name, fasta, thread, logFile, run = True):
    logFile.write(f"\n\n======================== Norm and Biallelic and Annotate ({name}) ========================\n\n")
    for i in range(begin, end + 1):
        vcfFileWithChr = vcfFile.replace("*", str(i))

        execute(f"{bcftools} norm -m -any {vcfFileWithChr} --threads {thread} -Oz -o {folder}/{name}_NormMinus_chr{i}.vcf.gz", logFile, run)
        execute(f"{bcftools} index {folder}/{name}_NormMinus_chr{i}.vcf.gz --threads {thread}", logFile, run)
        execute(f"{bcftools} norm -Oz --check-ref s -f {fasta} --threads {thread} -o {folder}/{name}_Norm_chr{i}.vcf.gz "
                f"{folder}/{name}_NormMinus_chr{i}.vcf.gz", logFile, run)
        execute(f"{bcftools} index {folder}/{name}_Norm_chr{i}.vcf.gz --threads {thread}", logFile, run)
        execute(f"{bcftools} view --types snps -m 2 -M 2 {folder}/{name}_Norm_chr{i}.vcf.gz -Oz "
                f"-o {folder}/{name}_Biallelic_chr{i}.vcf.gz --threads {thread}", logFile, run)
        execute(f"{bcftools} index {folder}/{name}_Biallelic_chr{i}.vcf.gz --threads {thread}", logFile, run)

        execute(f"{bcftools} annotate -x 'ID' -I +'%CHROM:%POS:%REF:%ALT' {folder}/{name}_Biallelic_chr{i}.vcf.gz -Oz "
                f"-o {folder}/{name}_Annotated_chr{i}.vcf.gz --threads {thread}", logFile, run)
        execute(f"{bcftools} index {folder}/{name}_Annotated_chr{i}.vcf.gz --threads {thread}", logFile, run)

    return f"{folder}/{name}_Annotated_chr*.vcf.gz"


def filterReference(reference, correspondenceList, bcftools, folder, begin, end, thread, logFile, run = True):
    logFile.write(f"\n\n============================= Filter Reference =============================\n\n")

    inputFile = open(correspondenceList)
    fileToFilter = open(f"{folder}/listToKeep.txt", "w")

    for line in inputFile:
        split = line.strip().split()
        fileToFilter.write(f"{split[0]}\n")

    fileToFilter.close()

    for i in range(begin, end + 1):
        referenceWithChr = reference.replace('*', str(i))
        execute(
            f"{bcftools} view {referenceWithChr} -S {folder}/listToKeep.txt -Oz -o {folder}/ReferenceFiltered_chr{i}.vcf.gz --threads {thread}",
            logFile, run)
        execute(f"{bcftools} index {folder}/ReferenceFiltered_chr{i}.vcf.gz --threads {thread}", logFile, run)

    return f"{folder}/ReferenceFiltered_chr*.vcf.gz"


def convertToVCF(plinkFile, plink, folder, name, begin, end, logFile, run = True):
    logFile.write(f"\n\n============================= VCF convertion =============================\n\n")
    for i in range(begin, end + 1):
        plinkWithChromosome = plinkFile.replace("*", str(i))
        execute(f'{plink} --bfile {plinkWithChromosome} --recode vcf-iid --out {folder}/{name}_chr{i} --output-chr chrMT', logFile, run)
    return f'{folder}/{name}_chr*.vcf'

def bgzip(vcf, bgzip, bcftools, begin, end, logFile, run = True):
    logFile.write(f"\n\n============================= VCF convertion =============================\n\n")
    for i in range(begin, end+1):
        vcfWithChromosome = vcf.replace("*", str(i))
        execute(f'{bgzip} {vcfWithChromosome}', logFile, run)
        execute(f'{bcftools} index {vcfWithChromosome}.gz', logFile, run)

    return f"{vcf}.gz"

def execute(line, logFile, run = True):
    print(f"{line}")
    if run:
        os.system(f"{line}")
    logFile.write(f"{line}\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='GWAS (and maybe XWAS) with Local Ancestry')

    requiredGeneral = parser.add_argument_group("Required arguments for all steps")
    requiredGeneral.add_argument('-i', '--input',
                                 help='Input file with the samples to infer the LA with chromosome replaced by \'*\'. '
                                      'If it is PLINK file (ex.bed, ex.bim, ex.fam), put just the prefix (-i ex). '
                                      'If it is VCF, put the whole name (-i ex.vcf or ex.vcf.gz)',
                                 required=True)
    requiredGeneral.add_argument('-I', '--inputImputed', help='Input file Imputed', required=True)
    requiredGeneral.add_argument('-o', '--outputName', help='Name of output name', required=True)
    requiredGeneral.add_argument('-O', '--outputFolder', help='Name of output folder', required=True)
    requiredGeneral.add_argument('-b', '--begin', help='First chromosome (default = 1)', default=1, type=int)
    requiredGeneral.add_argument('-e', '--end', help='Last chromosome  (default = 22)', default=22, type=int)

    requiredPhasing = parser.add_argument_group("Required arguments for Local Ancestry using RFMix1")
    requiredPhasing.add_argument('-l', '--referenceLA',
                                 help='Path to the reference database to perform the LA with chromosome replaced by *. '
                                      'If it is PLINK file (ex.bed, ex.bim, ex.fam), put just the prefix (-i ex). '
                                      'If it is VCF (ex.vcf or ex.vcf.gz), put the whole name (-i ex.vcf or ex.vcf.gz)',
                                 required=True)
    requiredPhasing.add_argument('-c', '--correspondence',
                                 help='Correspondence of reference database with ancestry. We will use just the '
                                      'samples in this file to perform the LA',
                                 required=True)
    requiredPhasing.add_argument('-p', '--referencePhase',
                                 help='Path to the reference database to phase the reference LA and Target with '
                                      'chromosome replaced by *. Eagle uses vcf files as reference panel',
                                 required=False, default='')
    requiredPhasing.add_argument('-s', '--setSize', help='Set size to run the RFMix v1 (default = 100)',
                                 required=True, default=100, type=int)
    requiredPhasing.add_argument('-g', '--geneticMap', help='Folder with genetic maps with chromosome replaced by *', required=True)
    requiredPhasing.add_argument('-f', '--fasta', help='Fasta of genome reference (to fix REF/ALT alleles)',
                                 required=True)

    programs = parser.add_argument_group("Programs")
    programs.add_argument('-E', '--eagle', help='Path to eagle', required=False, default="eagle")
    programs.add_argument('-B', '--bcftools', help='Path to bcftools', required=False, default="bcftools")
    programs.add_argument('-P', '--plink', help='Path to PLINK', required=False, default="plink")
    programs.add_argument('-R', '--rfmix1', help='Path to RFMix1', required=False, default="rfmix1")
    programs.add_argument('-Z', '--bgzip', help='Path to bgzip', required=False, default="bgzip")

    optional = parser.add_argument_group("Optional arguments")
    optional.add_argument('-t', '--threads', help='Number of threads (default = 20)', default=20, type=int,
                          required=False)

    args = parser.parse_args()

    print(f"Creating the output folder {args.outputFolder}")
    os.system(f"mkdir {args.outputFolder}")

    logFile = open(f"{args.outputFolder}/{args.outputName}.log", 'w')
    execute(f"mkdir {args.outputFolder}/Reference", logFile)
    execute(f"mkdir {args.outputFolder}/Target", logFile)

    # Aquecimento : mudando Bed bim fam para VCF
    if "vcf" not in args.input:
        target = convertToVCF(args.input, args.plink, f'{args.outputFolder}/Target', "Target", args.begin, args.end, logFile)
    else:
        target = args.target
    if "vcf" not in args.referenceLA:
        referenceLA = convertToVCF(args.referenceLA, args.plink, f'{args.outputFolder}/Reference', "Reference", args.begin, args.end, logFile)
    else:
        referenceLA = args.referenceLA

    # Se não VCF.gz, faça ser vcf.gz com index
    if "gz" not in target:
        target = bgzip(target, args.bgzip, args.bcftools, args.begin, args.end, logFile)
    if "gz" not in referenceLA:
        referenceLA = bgzip(referenceLA, args.bgzip, args.bcftools, args.begin, args.end, logFile)

    #Filter the reference
    referenceVCF = filterReference(referenceLA, args.correspondence, args.bcftools, f"{args.outputFolder}/Reference",
                                   args.begin, args.end, args.threads, logFile)


    #Merge process
    referenceVCF = normAndBiallelic(referenceVCF, f"{args.outputFolder}/Reference", args.bcftools, args.begin, args.end,
                                    "Reference", args.fasta, args.threads, logFile)

    targetVCF = normAndBiallelic(target, f"{args.outputFolder}/Target", args.bcftools, args.begin, args.end,
                                 "Target", args.fasta, args.threads, logFile)

    execute(f"mkdir {args.outputFolder}/Merged", logFile)
    allSamplesVCF = mergeVCFs(referenceVCF, targetVCF, args.bcftools, f'{args.outputFolder}/Merged', args.begin,
                              args.end, args.threads, logFile)

    execute(f"mkdir {args.outputFolder}/Phased", logFile)
    allSamplesPhased = phaseWithEagle(allSamplesVCF, args.referencePhase, f'{args.outputFolder}/Phased', args.begin,
                              args.end, args.threads, args.geneticMap, args.eagle, logFile)

    execute(f"mkdir {args.outputFolder}/PLINK", logFile)

    PCA = calculatePCAWithPlink(targetVCF, args.plink, f'{args.outputFolder}/PLINK', args.begin,
                                args.end, args.outputName, logFile)

    setsFiles = createSets(PCA, args.setSize, f'{args.outputFolder}/PLINK', args.outputName, logFile)