import argparse
import gzip
import time
import os

def splitGeneticMapByChromosome(geneticMap, folder, logFile):
    execute(f"mkdir {folder}/GeneticMaps/", logFile)
    if ".gz" in geneticMap:
        inputFile = gzip.open(geneticMap)
    else:
        inputFile = open(geneticMap)

    fileDict = {}

    header = True
    for line in inputFile:
        line = line.decode('utf-8')
        if header:
            headerLine = "pos chr cM\n"
            header = False
        else:
            split = line.strip().split()
            if split[0] not in fileDict:
                fileDict[split[0]] = open(f'{folder}/GeneticMaps/geneticMap_chr{split[0]}.txt', 'w')
                fileDict[split[0]].write(headerLine)
            fileDict[split[0]].write(f"{split[1]} {split[0]} {split[3]}\n")

    for chrom in fileDict:
        fileDict[chrom].close()

    return f'{folder}/GeneticMaps/geneticMap_chr*.txt'


def runRFMixSequentially(setVCF, originalVCF, numSet, begin, end, folder, name, threads, rfmix1, python2, python3, plink, VCF2RFMix,
                         RFMix1ToRFMix2, correspondence, geneticMap, setList, allelesRephased2VCF, doNotRephase, logFile, run = True):

    execute(f"mkdir {folder}/RFMix1_Outputs/", logFile, run)
    execute(f"mkdir {folder}/RFMix1_Inputs/", logFile, run)
    execute(f"mkdir {folder}/RFMix2/", logFile, run)

    for chrom in range(begin, end + 1):
        for setID in range(0, numSet+1):
            setVCFWithNumbers = setVCF.replace("#", str(setID)).replace("*", str(chrom))
            geneticMapWithChrom = geneticMap.replace("*", str(chrom))
            command = f"{python3} {VCF2RFMix} -v {setVCFWithNumbers} -c {correspondence} -C {chrom} " \
                      f"-o {folder}/RFMix1_Inputs/{name}_chrom{chrom}_set{setID} -m {geneticMapWithChrom} -p {plink}"
            execute(command, logFile, run)

            alleles = f'{folder}/RFMix1_Inputs/{name}_chrom{chrom}_set{setID}_alleles'
            classes = f'{folder}/RFMix1_Inputs/{name}_chrom{chrom}_set{setID}_classes'
            location = f'{folder}/RFMix1_Inputs/{name}_chrom{chrom}_set{setID}_location'

            type = "PopPhased"
            if doNotRephase:
                type = "TrioPhased"

            command = f"{python2} {rfmix1} {type} {alleles} {classes} {location} " \
                      f"-o {folder}/RFMix1_Outputs/Output_{name}_chrom{chrom}_set{setID} --num-threads {threads} -e 2 " \
                      f"-w 0.2 --forward-backward --skip-check-input-format --succinct-output"
            execute(command, logFile, run)

        originalVCFWithChrom = originalVCF.replace("*", str(chrom))
        FB = f'{folder}/RFMix1_Outputs/Output_{name}_chrom{chrom}_set\*.2.ForwardBackward.txt'
        SNPPerWindow = f'{folder}/RFMix1_Outputs/Output_{name}_chrom{chrom}_set0.2.SNPsPerWindow.txt'
        classes = f'{folder}/RFMix1_Inputs/{name}_chrom{chrom}_set\*_classes'

        command = f"{python3} {RFMix1ToRFMix2} -v {originalVCFWithChrom} -s {setList} " \
                  f"-m {correspondence} -b 0 -e {numSet} -C {classes} -F {FB} -c {chrom} -g {geneticMapWithChrom} " \
                  f"-o {folder}/RFMix2/{name}_{chrom} -S {SNPPerWindow}"
        execute(command, logFile, run)



        FB = f'{folder}/RFMix1_Outputs/Output_{name}_chrom{chrom}_set\*.2.ForwardBackward.txt'
        AR = f'{folder}/RFMix1_Outputs/Output_{name}_chrom{chrom}_set\*.allelesRephased2.txt'
        SNPPerWindow = f'{folder}/RFMix1_Outputs/Output_{name}_chrom{chrom}_set\*.2.SNPsPerWindow.txt'
        classes = f'{folder}/RFMix1_Inputs/{name}_chrom{chrom}_set\*_classes'
        execute(f"mkdir {folder}/VCFRephased", logFile)

        setVCFWithChrom = setVCF.replace("*", str(chrom)).replace("#", "\*")
        command = f"{python3} {allelesRephased2VCF} -v {setVCFWithChrom} -o {folder}/VCFRephased/{name}_{chrom}_Rephased.vcf " \
                  f"-c {correspondence} -s {setList} -b 0 -e {numSet} -S {SNPPerWindow} -C {classes} -F {FB} -A {AR}"
        execute(command, logFile, run)

def countFileLine(inputFile):
    file = open(inputFile)
    numLine = 0
    for line in file:
        numLine = numLine + 1
    return numLine

def readModelFile(modelFile):
    inputFile = open(modelFile)
    fileLines = ""

    for line in inputFile:
        fileLines = fileLines+line

    return fileLines

def generateListOfFlags(name, begin, end, numSet):
    listOfFlags = []
    for chrom in range(begin, end + 1):
        for setID in range(0, numSet+1):
            listOfFlags.append(f"{name}_{chrom}_{setID}.txt")

    return listOfFlags

def runRFMixBot(setVCF, originalVCF, numSet, begin, end, folder, name, rfmix1, python2, python3, plink, VCF2RFMix,
                RFMix1ToRFMix2, correspondence, geneticMap, setList, numJobs, queueCheck, queueSubmit, modelFile, memory,
                cores, allelesRephased2VCF, doNotRephase, logFile, run = True):


    execute(f"mkdir {folder}/RFMix1_Outputs/", logFile, run)
    execute(f"mkdir {folder}/RFMix1_Inputs/", logFile, run)
    execute(f"mkdir {folder}/ToSubmit/", logFile, run)
    execute (f"rm {folder}/ToSubmit/*", logFile, run)
    execute(f"mkdir {folder}/Flag/", logFile, run)
    execute(f"mkdir {folder}/RFMix2/", logFile, run)


    modelToUse = readModelFile(modelFile)
    print(modelToUse)

    allFlags = generateListOfFlags(name, begin, end, numSet)

    for chrom in range(begin, end + 1):
        for setID in range(0, numSet+1):
            putToRun = False
            while not putToRun:
                time.sleep(1)
                os.system(f'{queueCheck} > batch.txt')
                numLine = countFileLine('batch.txt')

                if numLine < numJobs:
                    modelToSave = modelToUse.replace("%%%NAME%%%", f"LA_{chrom}_{setID}").replace("%%%MEM%%%", memory).replace("%%%CORES%%%", cores)
                    print(modelToSave)
                    fileToSubmit = open(f"{folder}/ToSubmit/Submission_chrom{chrom}_set{setID}.sh", 'w')
                    fileToSubmit.write(f'{modelToSave}\n')

                    setVCFWithNumbers = setVCF.replace("#", str(setID)).replace("*", str(chrom))
                    geneticMapWithChrom = geneticMap.replace("*", str(chrom))
                    command = f"{python3} {VCF2RFMix} -v {setVCFWithNumbers} -c {correspondence} -C {chrom} " \
                              f"-o {folder}/RFMix1_Inputs/{name}_chrom{chrom}_set{setID} -m {geneticMapWithChrom} -p {plink}"

                    fileToSubmit.write(f'{command}\n')

                    alleles = f'{folder}/RFMix1_Inputs/{name}_chrom{chrom}_set{setID}_alleles'
                    classes = f'{folder}/RFMix1_Inputs/{name}_chrom{chrom}_set{setID}_classes'
                    location = f'{folder}/RFMix1_Inputs/{name}_chrom{chrom}_set{setID}_location'

                    type = "PopPhased"
                    if doNotRephase:
                        type = "TrioPhased"

                    command = f"{python2} {rfmix1} {type} {alleles} {classes} {location} " \
                              f"-o {folder}/RFMix1_Outputs/Output_{name}_chrom{chrom}_set{setID} --num-threads {cores} -e 2 " \
                              f"-w 0.2 --forward-backward --skip-check-input-format --succinct-output"
                    fileToSubmit.write(f'{command}\n')
                    fileToSubmit.write(f'> {folder}/Flag/{name}_{chrom}_{setID}.txt\n')
                    fileToSubmit.close()
                    command = f"{queueSubmit} {folder}/ToSubmit/Submission_chrom{chrom}_set{setID}.sh"

                    putToRun = True
                    execute(command, logFile, run)

    allDone = False
    while not allDone:
        time.sleep(5)
        folderFiles = os.listdir(f"{folder}/Flag/")
        allDone = True
        for flag in allFlags:
            if flag not in folderFiles:
                allDone = False

    for chrom in range(begin, end + 1):
        geneticMapWithChrom = geneticMap.replace("*", str(chrom))
        originalVCFWithChrom = originalVCF.replace("*", str(chrom))
        FB = f'{folder}/RFMix1_Outputs/Output_{name}_chrom{chrom}_set\*.2.ForwardBackward.txt'
        SNPPerWindow = f'{folder}/RFMix1_Outputs/Output_{name}_chrom{chrom}_set0.2.SNPsPerWindow.txt'
        classes = f'{folder}/RFMix1_Inputs/{name}_chrom{chrom}_set\*_classes'

        command = f"{python3} {RFMix1ToRFMix2} -v {originalVCFWithChrom} -s {setList} " \
                  f"-m {correspondence} -b 0 -e {numSet} -C {classes} -F {FB} -c {chrom} -g {geneticMapWithChrom} " \
                  f"-o {folder}/RFMix2/{name}_{chrom} -S {SNPPerWindow}"
        execute(command, logFile, run)

        FB = f'{folder}/RFMix1_Outputs/Output_{name}_chrom{chrom}_set\*.2.ForwardBackward.txt'
        AR = f'{folder}/RFMix1_Outputs/Output_{name}_chrom{chrom}_set\*.allelesRephased2.txt'
        SNPPerWindow = f'{folder}/RFMix1_Outputs/Output_{name}_chrom{chrom}_set\*.2.SNPsPerWindow.txt'
        classes = f'{folder}/RFMix1_Inputs/{name}_chrom{chrom}_set\*_classes'
        execute(f"mkdir {folder}/VCFRephased", logFile)

        setVCFWithChrom = setVCF.replace("*", str(chrom)).replace("#", "\*")
        command = f"{python3} {allelesRephased2VCF} -v {setVCFWithChrom} -o {folder}/VCFRephased/{name}_{chrom}_Rephased.vcf " \
                  f"-c {correspondence} -s {setList} -b 0 -e {numSet} -S {SNPPerWindow} -C {classes} -F {FB} -A {AR}"
        execute(command, logFile, run)


def separateSets(VCF, clusteringFile, correspondence, bcftools, folder, name, begin, end, threads, logFile, run = True):
    inputFile = open(correspondence)
    dictRef = {}
    dictSet = {}

    for line in inputFile:
        ID, anc = line.strip().split()
        dictRef[ID] = anc
    inputFile.close()

    header = True
    inputFile = open(clusteringFile)
    for line in inputFile:
        if header:
            header = False
        else:
            split = line.strip().split()

            if split[-1] not in dictSet:
                dictSet[split[-1]] = []
            dictSet[split[-1]].append(split[0])
    inputFile.close()

    maxSet = 0
    for setSamples in dictSet:
        if int(setSamples) > maxSet:
            maxSet = int(setSamples)
        fileToExtract = open(f"{folder}/ExtractSet{setSamples}.txt", 'w')

        for ind in dictSet[setSamples]:
            fileToExtract.write(f"{ind}\n")
        for ind in dictRef:
            fileToExtract.write(f"{ind}\n")
        fileToExtract.close()
        for chrom in range(begin, end + 1):
            filterVCF(VCF, f"{folder}/ExtractSet{setSamples}.txt", bcftools, folder, f"{name}_set{setSamples}", chrom,
                      threads, logFile, run)

    return f"{folder}/{name}_set#_chr*.vcf.gz", maxSet, f"{folder}/ExtractSet\*.txt"




def createSets(eigenvec, setSize, outputFolder, outputName, numPCs, logFile):

    #Code from https://stackoverflow.com/questions/5452576/k-means-algorithm-variation-with-equal-cluster-size?rq=1
    from sklearn.cluster import KMeans
    from scipy.spatial.distance import cdist
    from scipy.optimize import linear_sum_assignment
    import numpy as np
    import pandas as pd

    namesList = ['IID', 'FID']

    for i in range(1, int(numPCs)+1):
        namesList.append(f'PC{i}')

    X = pd.read_table(eigenvec, header=None, names=namesList, sep=" ", index_col="IID", usecols=['IID', 'PC1', 'PC2'])

    n_clusters = int(np.ceil(len(X)/setSize))
    kmeans = KMeans(n_clusters)
    kmeans.fit(X)
    centers = kmeans.cluster_centers_
    centers = centers.reshape(-1, 1, X.shape[-1]).repeat(setSize, 1).reshape(-1, X.shape[-1])
    distance_matrix = cdist(X, centers)
    clusters = linear_sum_assignment(distance_matrix)[1]//setSize

    X['ClusterID'] = clusters
    X.to_csv(f"{outputFolder}/{outputName}_PCA_Cluster.csv", sep = '\t')

    return f"{outputFolder}/{outputName}_PCA_Cluster.csv"

def phaseWithEagle(allSamplesVCF, referencePhase, outputFolder, begin, end, threads, geneticMap, eagle, bcftools, logFile, run = True):
    for i in range(begin, end + 1):
        samplesWithChr = allSamplesVCF.replace('*', str(i))
        vcfRefWithChr = referencePhase.replace('*', str(i))
        command = f'{eagle} --vcfTarget {samplesWithChr} --outPrefix {outputFolder}/Merged_Phased_{i} --numThreads {threads}' \
                  f' --geneticMapFile {geneticMap} --keepMissingPloidyX'
        if referencePhase != '':
            command = f'{command} --vcfRef {vcfRefWithChr} --allowRefAltSwap'
        execute(command, logFile, run)
        execute(f"{bcftools} index {outputFolder}/Merged_Phased_{i}.vcf.gz --threads {threads}", logFile, run)
    return f'{outputFolder}/Merged_Phased_*.vcf.gz'


def calculatePCAWithPlink(vcf, plink, folder, begin, end, outputName, numPCs, logFile, run = True):
    logFile.write(f"\n\n============================= PLINK PCA =============================\n\n")

    fileMerge = open(f'{folder}/ListToConcat.txt', 'w')
    for i in range(begin, end + 1):
        vcfWithChr = vcf.replace("*", str(i))
        execute(f"{plink} --vcf {vcfWithChr} --double-id --make-bed --out {folder}/{outputName}_chr{i}", logFile, run)
        fileMerge.write(f"{folder}/{outputName}_chr{i}\n")
    fileMerge.close()

    execute(f"{plink} --merge-list {folder}/ListToConcat.txt --make-bed --out {folder}/{outputName}_All", logFile, run)
    execute(f"{plink} --bfile {folder}/{outputName}_All --indep-pairwise 200 50 0.2 --out {folder}/{outputName}_LD",
            logFile, run)
    execute(f"{plink} --bfile {folder}/{outputName}_All --extract {folder}/{outputName}_LD.prune.in --make-bed --out "
            f"{folder}/{outputName}_LD_Prunned", logFile, run)
    execute(f"{plink} --bfile {folder}/{outputName}_LD_Prunned --pca {numPCs} --out {folder}/{outputName}_PCA", logFile, run)

    return f'{folder}/{outputName}_PCA.eigenvec'


def mergeVCFs(vcf1, vcf2, bcftools, folder, begin, end, thread, logFile, run = True):
    logFile.write(f"\n\n============================= Isec and Merge =============================\n\n")
    for i in range(begin, end + 1):
        vcf1WithChr = vcf1.replace("*", str(i))
        vcf2WithChr = vcf2.replace("*", str(i))

        execute(f"{bcftools} isec -n=2 {vcf1WithChr} {vcf2WithChr} -o {folder}/CommonVariants_chr{i} --threads {thread}", logFile, run)
        execute(f"{bcftools} merge -R {folder}/CommonVariants_chr{i} {vcf1WithChr} {vcf2WithChr} -Oz "
                f"-o {folder}/Merged_chr{i}.vcf.gz", logFile, run)

        #execute(f"{bcftools} merge -R {folder}/CommonVariants_chr{i} {vcf1WithChr} {vcf2WithChr} -Oz "
        #        f"-o {folder}/Merged_chr{i}.vcf.gz --threads {thread}", logFile, run)
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


def filterVCF(VCF, fileToFilter, bcftools, folder, name, chrom, thread, logFile, run):
    VCFWithChr = VCF.replace('*', str(chrom))
    execute(f"{bcftools} view {VCFWithChr} -S {fileToFilter} -Oz -o {folder}/{name}_chr{chrom}.vcf.gz --threads {thread}",
        logFile, run)
    execute(f"{bcftools} index {folder}/{name}_chr{chrom}.vcf.gz --threads {thread}", logFile, run)


def filterReference(reference, correspondenceList, bcftools, folder, begin, end, thread, logFile, run = True):
    logFile.write(f"\n\n============================= Filter Reference =============================\n\n")

    inputFile = open(correspondenceList)
    fileToFilter = open(f"{folder}/listToKeep.txt", "w")

    for line in inputFile:
        split = line.strip().split()
        fileToFilter.write(f"{split[0]}\n")

    fileToFilter.close()

    for i in range(begin, end + 1):
        filterVCF(reference, f"{folder}/listToKeep.txt", bcftools, folder, "ReferenceFiltered", i, thread, logFile, run)

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
    parser = argparse.ArgumentParser(description='Run Local Ancestry with RFMix v1')

    requiredGeneral = parser.add_argument_group("Required arguments for all steps")
    requiredGeneral.add_argument('-i', '--input',
                                 help='Input file with the samples to infer the LA with chromosome replaced by \'*\'. '
                                      'If it is PLINK file (ex.bed, ex.bim, ex.fam), put just the prefix (-i ex). '
                                      'If it is VCF, put the whole name (-i ex.vcf or ex.vcf.gz)',
                                 required=True)
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
    requiredPhasing.add_argument('-g', '--geneticMap', help='Folder with genetic map (downloaded from Eagle website, all '
                                                            'chromosome in the same file)', required=True)
    requiredPhasing.add_argument('-f', '--fasta', help='Fasta of genome reference (to fix REF/ALT alleles)',
                                 required=True)


    requiredRegression = parser.add_argument_group("Required arguments for Regression")
    requiredRegression.add_argument('-n', '--numPCs', help='Number of PCs to use in the regression',
                                 required=True)

    programs = parser.add_argument_group("Programs")
    programs.add_argument('-E', '--eagle', help='Path to eagle', required=False, default="eagle")
    programs.add_argument('-B', '--bcftools', help='Path to bcftools', required=False, default="bcftools")
    programs.add_argument('-P', '--plink', help='Path to PLINK', required=False, default="plink")
    programs.add_argument('-R', '--rfmix1', help='Path to RFMix1', required=False, default="rfmix1")
    programs.add_argument('-Z', '--bgzip', help='Path to bgzip', required=False, default="bgzip")
    programs.add_argument('-Y', '--python3', help='Path Python3 interpreter', required=False, default="python3")
    programs.add_argument('-y', '--python2', help='Path Python2 interpreter', required=False, default="python2")
    programs.add_argument('-V', '--VCF2RFMix', help='Path VCF2RFMix script', required=False, default="VCF2RFMix.py")
    programs.add_argument('-r', '--RFMix1ToRFMix2', help='Path RFMix1ToRFMix2 script', required=False, default="RFMix1ToRFMix2.py")
    programs.add_argument('-A', '--allelesRephased2VCF', help='Path allelesRephased2VCF script', required=False, default="allelesRephased2VCF.py")

    optional = parser.add_argument_group("Optional arguments")
    optional.add_argument('-t', '--threads', help='Number of threads (default = 20)', default=20, type=int,
                          required=False)
    optional.add_argument('-d', '--doNotRephase', help='Do not rephase using RFMix v1', default=False, action="store_true")

    botArgs = parser.add_argument_group("Bot Mode", description= "This mode was implemented to keep a minimun number of "
                                                                 "jobs in the submission queue. We implemented this "
                                                                 "feature in a slurm HPC. The goal is look the queue and "
                                                                 "automatically add a job if the number of jobs is lower "
                                                                 "than the requested. If not activated, our pipeline will"
                                                                 "proceed the RFMix inference sequentially.")
    botArgs.add_argument('-j', '--jobs', help='Activate the bot mode and set the number of jobs on the queue',
                         required=False, default= -1, type=int)
    botArgs.add_argument('-q', '--queueCheck', help='Queue check command', required=False, default="squeue")
    botArgs.add_argument('-Q', '--queueSubmit', help='Command to submit a job on the queue', required=False, default="sbatch")
    botArgs.add_argument('-m', '--model', help='Model of the submission file to be changed.', required=False)
    botArgs.add_argument('-M', '--memory', help='Memory to be requested in the queue submission', required=False)
    botArgs.add_argument('-C', '--cores', help='Number of cores to be requested in the queue submission', required=False)


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
        target = args.input
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

    #Phasing

    if not args.doNotRephase:
        execute(f"mkdir {args.outputFolder}/Phased", logFile)
        allSamplesPhased = phaseWithEagle(allSamplesVCF, args.referencePhase, f'{args.outputFolder}/Phased', args.begin,
                                  args.end, args.threads, args.geneticMap, args.eagle, args.bcftools, logFile)
    else:
        allSamplesPhased = allSamplesVCF

    #Sets
    execute(f"mkdir {args.outputFolder}/PLINK", logFile)
    PCA = calculatePCAWithPlink(targetVCF, args.plink, f'{args.outputFolder}/PLINK', args.begin,
                                args.end, args.outputName, args.numPCs, logFile)

    clusteringFile = createSets(PCA, args.setSize, f'{args.outputFolder}/PLINK', args.outputName, args.numPCs, logFile)

    execute(f"mkdir {args.outputFolder}/Sets", logFile)
    setVCF, numSet, setList = separateSets(allSamplesPhased, clusteringFile, args.correspondence, args.bcftools,
                                  f"{args.outputFolder}/Sets", args.outputName, args.begin, args.end, args.threads, logFile)

    execute(f"mkdir {args.outputFolder}/RFMix", logFile)
    geneticMapSplit = splitGeneticMapByChromosome(args.geneticMap, f"{args.outputFolder}/RFMix", logFile)
    print(f"I have {args.jobs}")
    if args.jobs <= 0:
        runRFMixSequentially(setVCF, allSamplesPhased, numSet, args.begin, args.end, f"{args.outputFolder}/RFMix",
                             args.outputName, args.threads, args.rfmix1, args.python2, args.python3, args.plink,
                             args.VCF2RFMix, args.RFMix1ToRFMix2, args.correspondence, geneticMapSplit, setList,
                             args.allelesRephased2VCF, args.doNotRephase, logFile)
    else:
        runRFMixBot(setVCF, allSamplesPhased, numSet, args.begin, args.end, f"{args.outputFolder}/RFMix", args.outputName,
                    args.rfmix1, args.python2, args.python3, args.plink, args.VCF2RFMix, args.RFMix1ToRFMix2,
                    args.correspondence, geneticMapSplit, setList, args.jobs, args.queueCheck, args.queueSubmit.replace("\"", ""),
                    args.model, args.memory, args.cores, args.allelesRephased2VCF, args.doNotRephase, logFile)


