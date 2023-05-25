import os
import sys
import gzip
import argparse
version = "alpha"

def basicInfoVCF(ancIDPerSetDict, refID):
    string = "##fileformat=VCFv4.1\n"
    string = string + '##FILTER=<ID=PASS,Description="All filters passed\">\"\n'
    string = string + '##Ancestry IDs:'
    for id in ancIDPerSetDict[refID]:
        string = string + f' {id}:{ancIDPerSetDict[refID][id]}'
    string = string+ "\n"

    return string

def metaInfoVCF(rawLine):
    string = f"##Alleles_Rephase_version={version}\n"
    string = string + f"##Alleles_Rephase_Command_Line={rawLine}\n"
    return string

def infoVCF(ancIDPerSetDict, refID):
    string = f'##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency for alternative allele">\n'
    string = string + f'##INFO=<ID=COUNT,Number=A,Type=Integer,Description="Allele count for alternative allele">\n'
    for anc in ancIDPerSetDict[refID]:
        string = string + f'##INFO=<ID=COUNT_{anc},Number=A,Type=Integer,Description="Allele count alternative allele with ID {anc} ({ancIDPerSetDict[refID][anc]}) ancestry">\n'
        string = string + f'##INFO=<ID=AF_{anc},Number=A,Type=Float,Description="Allele frequency alternative allele with ID {anc} ({ancIDPerSetDict[refID][anc]}) ancestry">\n'

    return string

def infoFormatVCF():
    string = f'##FORMAT=<ID=GT,Number=1,Type=String,Description="Rephased genotype">\n'
    string = string + f'##FORMAT=<ID=FA,Number=1,Type=Integer,Description="First allele ancestry">\n'
    string = string + f'##FORMAT=<ID=SA,Number=1,Type=Integer,Description="Second allele ancestry">\n'
    string = string + f'##FORMAT=<ID=FAPP,Number=1,Type=String,Description="First allele posteriori probability">\n'
    string = string + f'##FORMAT=<ID=SAPP,Number=1,Type=String,Description="Second allele posteriori probability">\n'

    return string

def headerVCF(setDict):
    string = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    for set in setDict:
        for ind in setDict[set]:
            string = string + f"\t{ind}"

    return string + "\n"

def readSetFile(set, begin, end, references):
    dictSetInd = {}
    for num in range(begin, end + 1):
        dictSetInd[num] = []
        setWithChr = set.replace("*", str(num))

        setFile = open(setWithChr)
        for line in setFile:
            line = line.strip()
            if line not in references:
                dictSetInd[num].append(line)

    return dictSetInd

def getAncIDPerSet(vcf, setDict, begin, end, classes, referencesDict):
    setAncIDPerSet = {}
    for setNum in range(begin, end+1):
        setAncIDPerSet[setNum] = {}
        vcfWithSetNum = vcf.replace("*", str(setNum))
        classesWithSetNum = classes.replace("*", str(setNum))


        classOpen = open(classesWithSetNum)
        for line in classOpen:
            classFile = line.strip().split()
        classOpen.close()

        if ".gz" in vcfWithSetNum:
            inputFile = gzip.open(vcfWithSetNum)
            decode = True
        else:
            inputFile = open(vcfWithSetNum)
            decode = False

        headerFound = False
        for line in inputFile:
            if decode:
                line = line.decode('utf-8')
            if not headerFound:
                if "#CHROM" in line:
                    headerFound = True
                    header = line.strip().split()
                    for i in range(9, len(header)):
                        id = header[i]
                        indexClass = (i-9)*2
                        classInd = classFile[indexClass] #TODO DO X chromosome
                        if classInd != '0' and id not in setDict[setNum]:
                            anc = referencesDict[id]
                            if classInd not in setAncIDPerSet[setNum]:
                                setAncIDPerSet[setNum][classInd] = anc
                            elif anc != setAncIDPerSet[setNum][classInd]:
                                input("ERRORRRRRRRRRRRRRRRRRRRRR")
            else:
                break

    return setAncIDPerSet

#Check if all SNP per window has the same information
#If yes, return a list with the SNP per window
#If not die
def readSNPPerWindowFile(snpPerWindow, begin, end):
    refFile = open(snpPerWindow.replace("*", str(begin)))

    SPWList = []
    for i in range(begin, end+1):
        SPWList.append(open(snpPerWindow.replace("*", str(i))))

    listSnpPerWindow = []
    for line in refFile:
        numRef = int(line.strip())
        listSnpPerWindow.append(numRef)
        for i in range(begin+1, end+1):
            numAlt = int(SPWList[i].readline().strip())
            if numRef != numAlt:
                sys.exit("The SNP per window script is different for each set")
    return listSnpPerWindow

def splitLine(line, decode):
    if decode:
        return line.decode("utf-8").split()
    return line.split()

#Check if all VCF has the same variants
#If yes, return true
#If not, return false
def VCFHasTheSameVariants(vcf, begin, end):
    inputVCFList = []
    for i in range(begin, end+1):
        inputVCFName = vcf.replace("*", "1")
        if ".gz" in inputVCFName:
            inputVCFList.append(gzip.open(inputVCFName))
            decode = True
        else:
            inputVCFList.append(open(inputVCFName))
            decode = False

    #Jump the header
    for i in range(begin, end+1):
        lineCount = 1
        if decode:
            while "#CHROM" not in inputVCFList[i].readline().decode("utf-8"):
                lineCount = lineCount+1

        else:
            while "#CHROM" not in inputVCFList[i].readline():
                lineCount = lineCount+1

    sameInformation = True
    for line in inputVCFList[0]:
        splitRef = splitLine(line, decode)
        for i in range(begin+1, end+1):
            splitOther = splitLine(inputVCFList[i].readline(), decode)

            for j in range(0, 4):
                if splitRef[j] != splitOther[j]:
                    print(f"Diffence between set 0 and set {j} : {splitRef} vs {splitOther} -> index {j}")
                    return False

    for i in range(begin, end+1):
        inputVCFList[i].close()
    print(f"All VCFs has the same information")

    return True


def getInfoFromFBLine(baseFB, ancCount, FBSplit):
    APP = f"{FBSplit[baseFB]}"

    maxProb = FBSplit[baseFB]
    maxIndex = 1
    for j in range(baseFB+1, baseFB+ancCount):
        APP = APP + f",{FBSplit[j]}"
        if float(maxProb) < float(FBSplit[j]):
            maxProb = FBSplit[j]
            maxIndex = (j-baseFB)+1

    return maxIndex, APP




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='RFMix1 to RFMix2')

    required = parser.add_argument_group("Required arguments")
    required.add_argument('-v', '--vcf', help='VCF separated by set with chr number replaced by * ', required=True)  # VCF

    # IO
    required.add_argument('-o', '--output', help='Output prefix', required=True)  # out and VCF

    # Set
    required.add_argument('-c', '--correspondence', help='Correspondence between parental and IDs', required=True)
    required.add_argument('-s', '--set', help='Sets file with set number replaced by *', required=True)  # out and VCF
    required.add_argument('-b', '--begin', help='First set', required=True, type=int)  # out
    required.add_argument('-e', '--end', help='Last set', required=True, type=int)  # out

    # RFMix files
    required.add_argument('-C', '--classes', help='Classes files with set number replaced by *',
                          required=True)
    required.add_argument('-S', '--snpPerWindow', help='SNPsPerWindow RFMix1 with the set number replaced by *',
                          required=True)
    required.add_argument('-A', '--allelesRephased',
                          help='allelesRephased from RFMix1 with the set number replaced by *', required=True)
    required.add_argument('-F', '--forwardBackward',
                          help='forwardBackward from RFMix1 with the set number replaced by *', required=True)

    # Plink
    required = parser.add_argument_group("optional arguments")
    required.add_argument('-p', '--plink', help='Path to PLINK (default plink)', required=False, default='plink')
    required.add_argument('-X', '--XList',
                          help='List of individuals created by VCF2RFMix with flag -X or --Xmen with the set number replaced by *',
                          required=False, default='')

    args = parser.parse_args()
    rawLine = ' '.join(f'{k}={v}' for k, v in vars(args).items())

    correspondence = open(f'{args.correspondence}')
    referencesDict = {}

    for line in correspondence:
        split = line.strip().split()
        referencesDict[split[0]] = split[1]

    setDict = readSetFile(args.set, args.begin, args.end, referencesDict)
    ancIDPerSetDict = getAncIDPerSet(args.vcf, setDict, args.begin, args.end, args.classes, referencesDict)

    snpPerWindow = readSNPPerWindowFile(args.snpPerWindow, args.begin, args.end)

    setRef = False
    allSetsHaveTheSameReference = True
    for setID in ancIDPerSetDict:
        if not setRef:
            setRef = True
            refID = setID
        else:
            ancCount = 0
            for ID in ancIDPerSetDict[setID]:
                ancCount = ancCount + 1
                if ancIDPerSetDict[setID][ID] != ancIDPerSetDict[refID][ID]:
                    print(f"Diff ID ({ID}) between set Ref (ID: {refID}) -> {ancIDPerSetDict[refID][ID]} and set ID {setID} -> {ancIDPerSetDict[setID][ID]}")


    print(f"We have {ancCount} parental populations")

    if VCFHasTheSameVariants(args.vcf, args.begin, args.end):
        header = True
        outputVCF = open(args.output, 'w')

        #Open a VCF to get CHROM, POS, Alleles information
        inputVCFName = args.vcf.replace("*", str(args.begin))
        if ".gz" in inputVCFName:
            vcfFile = gzip.open(inputVCFName)
            decode = True
        else:
            vcfFile = open(inputVCFName)
            decode = False

        allelesRephasedList = []
        forwardBackwardList = []
        for i in range(args.begin, args.end+1):
            allelesRephasedList.append(open(args.allelesRephased.replace("*", str(i))))
            forwardBackwardList.append(open(args.forwardBackward.replace("*", str(i))))

        print("Output the VCF:")
        for line in vcfFile:
            if decode:
                line = line.decode("utf-8")
            if header:
                if "#CHROM" not in line:
                    pass
                else:
                    header = False
                    print("\tCreating header")
                    outputVCF.write(basicInfoVCF(ancIDPerSetDict, refID))
                    outputVCF.write(metaInfoVCF(rawLine))
                    outputVCF.write(infoVCF(ancIDPerSetDict, refID))
                    outputVCF.write(infoFormatVCF())
                    outputVCF.write(headerVCF(setDict))
                    countVariants = 0
                    totalVariants = 0
            else:
                if countVariants == 0:
                    countVariants = snpPerWindow.pop(0)
                    #print(f"This window has {countVariants} SNPs")
                    FBSplit = forwardBackwardList[args.begin].readline().strip().split()
                    for i in range(args.begin+1, args.end+1):
                        FBSplit = FBSplit + forwardBackwardList[i].readline().strip().split()


                split = line.strip().split()
                outputVCF.write(f"{split[0]}\t{split[1]}\t{split[2]}\t{split[3]}\t{split[4]}\t{split[5]}\t{split[6]}\t")


                allRephased = ""
                for i in range(args.begin, args.end+1):
                    allRephased = allRephased+allelesRephasedList[i].readline().strip()

                alleleDict = {}
                alleleDict["All"] = 0
                alleleDict["Alt"] = 0

                for i in range(1, ancCount+1):
                    alleleDict[i] = {}
                    alleleDict[i]["Alt"] = 0
                    alleleDict[i]["All"] = 0

                sampleInfo = ""
                for i in range(int(len(allRephased)/2)):
                    baseRephased = i*2
                    baseFB = 2*i*ancCount

                    GT = f"{allRephased[baseRephased]}|{allRephased[baseRephased+1]}"

                    # Get First ancestry and Posteriori probabilities related to them
                    FA, FAPP = getInfoFromFBLine(baseFB, ancCount, FBSplit)
                    # Count Ref and All alleles for first ancestry
                    if str(allRephased[baseRephased]) == "1":
                        alleleDict["Alt"] = alleleDict["Alt"]+1
                        alleleDict[FA]["Alt"] = alleleDict[FA]["Alt"] + 1
                    alleleDict[FA]["All"] = alleleDict[FA]["All"] + 1

                    # Get Second ancestry and Posteriori probabilities related to them
                    SA, SAPP = getInfoFromFBLine(baseFB + ancCount, ancCount, FBSplit)
                    # Count Ref and All alleles for second ancestry
                    if str(allRephased[baseRephased+1]) == "1":
                        alleleDict["Alt"] = alleleDict["Alt"]+1
                        alleleDict[SA]["Alt"] = alleleDict[SA]["Alt"] + 1
                    alleleDict[SA]["All"] = alleleDict[SA]["All"] + 1

                    alleleDict["All"] = alleleDict["All"] + 2

                    sampleInfo = sampleInfo + f"\t{GT}:{FA}:{SA}:{FAPP}:{SAPP}"

                AF = "{:.3f}".format(alleleDict["Alt"]/alleleDict["All"])
                infoLine= f"AF={AF};COUNT={alleleDict['Alt']}"
                for i in range(1, ancCount+1):
                    if alleleDict["All"] != 0:
                        AF = "{:.3f}".format(alleleDict[i]["Alt"]/alleleDict[i]["All"])
                        infoLine = infoLine + f";COUNT_{i}={alleleDict[i]['Alt']};AF_{i}={AF}"
                    else:
                        infoLine = infoLine + f";COUNT_{i}=NA;AF_{i}=NA"

                outputVCF.write(f"{infoLine}\tGT:FA:SA:FAPP:SAPP{sampleInfo}\n")
                countVariants = countVariants - 1
                totalVariants = totalVariants + 1

    print(f"Your VCF has {totalVariants} ({countVariants})")
    outputVCF.close()