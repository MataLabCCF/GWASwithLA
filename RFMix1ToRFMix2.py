import os
import sys
import vcfpy
import argparse

def getchar():
    c = sys.stdin.read(1)

def createMappingDict(mapping):
    file = open(mapping)

    dict = {}

    for line in file:
        split = line.strip().split()
        dict[split[0]] = split[1]

    return dict

def createClassesPerSetDict(set, classes, mappingDict, first, last, XList):
    print('Function: createClassesPerSetDict')
    indPerSet = {}
    classIDPerSet = {}
    pops = []
    XDict = {}
    for setNumber in range(first, last+1):
        if XList:
            XListChr = XList.replace('*', str(setNumber))
            XFile = open(XListChr)
            for line in XFile:
                split = line.strip().split()
                XDict[split[0]] = split[1]

        setWithChr = set.replace('*', str(setNumber))
        classesWithChr = classes.replace('*', str(setNumber))

        print(f'\tOpen {setWithChr} to create a dict with individual by set')
        setOpen = open(setWithChr)
        indList = []

        indPerSet[setNumber] = []

        for line in setOpen:
            cleanLine = line.strip()
            if cleanLine not in mappingDict:
                indPerSet[setNumber].append(cleanLine)
            indList.append(cleanLine)
            indList.append(cleanLine)

        classIDPerSet[setNumber] = {}

        print(f'\tOpen {classesWithChr} to create a dict with ID and anc by set')
        classesOpen = open(classesWithChr)
        classesList = classesOpen.readline().split()

        for i in range(0, len(classesList)):
            id = classesList[i]
            if id != '0':
                ind = indList[i]
                anc = mappingDict[ind]
                classIDPerSet[setNumber][anc] = id
                if anc not in pops:
                    pops.append(anc)

        if not XList:
            for ind in indList:
                XDict[ind] = "D"


    print('\tAll dicts were created')

    # Dict: ClassIDPerSet
    # Key: SetNumber
    #   Key: Population ID
    #       Value: ID (1,2,3,...)
    # ---------------------------
    # Dict: indPerSet
    # Key: SetNumber
    #   Value: list with all individuals in this set
    # ---------------------------
    # List: pops
    # List with all parental populations
    # ---------------------------

    return classIDPerSet, indPerSet, pops, XDict
def createHeader(popList):
    header = '#reference_panel_population:'
    for parental in popList:
        header = f'{header}\t{parental}'
    header = header + '\nchromosome\tphysical_position\tgenetic_position\tgenetic_marker_index'
    return header
def addToHeader(indPerSet, popList, setNum, header):
    for ind in indPerSet[setNum]:
        for i in range(1,3):
            for parental in popList:
                header = header+f'\t{ind}:::hap{i}:::{parental}'
    return header

def createHeaderWithMSP(popList):
    header = '#reference_panel_population:'
    for parental in popList:
        header = f'{header}\t{parental}'
    header = header + '\nchromosome\tphysical_position\tgenetic_position\tgenetic_marker_index'
    
    headerMSP = f'#Subpopulation order/codes:'
    for i in range(0, len(popList)):
        headerMSP = f'{headerMSP}{popList[i]}={i}\t'
    headerMSP = headerMSP + f'\n#chm\tspos\tepos\tsgpos\tegpos\tn snps'
    
    return header, headerMSP

def addToHeaderWithMSP(indPerSet, popList, setNum, header, headerMSP, XDict):
    for ind in indPerSet[setNum]:
        if XDict[ind] == "D":
            for i in range(1,3):
                for parental in popList:
                    header = header+f'\t{ind}:::hap{i}:::{parental}'
            headerMSP = headerMSP+f'\t{ind}\t{ind}.1'
        else:
            for i in range(1,2):
                for parental in popList:
                    header = header+f'\t{ind}:::hap{i}:::{parental}'
            headerMSP = headerMSP+f'\t{ind}'
    return header, headerMSP

def convertOutputRFMixWindow(vcf, mapping, fowardBackward, output, setPrefix, classes, first, last, chrom, plink, geneticMap, Xlist, snpPerWindow):
    mappingDict = createMappingDict(mapping)
    classesPerSet, indPerSet, popList, XDict = createClassesPerSetDict(setPrefix, classes, mappingDict, first, last, Xlist)

    print(f'Open the fowardBackward files to create a {output}.fb.tsv file')

    header, headerMSP = createHeaderWithMSP(popList)
    firstFile = True
    snpPerWindowReaded = False

    for setNum in range(first, last+1):
        header, headerMSP = addToHeaderWithMSP(indPerSet, popList, setNum, header,headerMSP, XDict)
        fowardBackwardWithSet = fowardBackward.replace('*', str(setNum))

        if not snpPerWindowReaded:
            SNPsPerWindowDict = {}
            snpPerWindowReaded = True
            nameSNPsPerWindow = snpPerWindow.replace('*', str(setNum))
            
            print(f'\tOpen the {nameSNPsPerWindow} file and store this data on a dict')
            
            fileSNPsPerWindow = open(nameSNPsPerWindow)
            
            window = 0
            for line in fileSNPsPerWindow:
                SNPsPerWindowDict[window] = int(line.strip())    
                window = window+1
            
        if firstFile:
            firstFile = False
            print(f'\tOpen the {fowardBackwardWithSet} file and create the first tempFile {output}_set{setNum}.fb.tsv')
            fileFB = open(fowardBackwardWithSet)
            fbOut = open(f'{output}_set{setNum}.fb.tsv', 'w')
            for line in fileFB:
                posterioris = line.strip().split()
                for i in range(0,len(posterioris), len(popList)):
                    for parental in popList:
                        classID = int(classesPerSet[setNum][parental]) - 1 #starts on 0
                        post = posterioris[classID+i]
                        fbOut.write(f'\t{post}')
                fbOut.write('\n')
            fileFB.close()
            fbOut.close()
            print(f'\tDone to file {fowardBackwardWithSet}')
        else:
            previous = setNum-1
            print(f'\tOpen the {fowardBackwardWithSet} file and create tempFile {output}_set{setNum}.fb.tsv')
            fileFB = open(fowardBackwardWithSet)
            fbPrevious = open(f'{output}_set{previous}.fb.tsv')
            fbOut = open(f'{output}_set{setNum}.fb.tsv', 'w')
            for line in fileFB:
                linePrevious = fbPrevious.readline().strip()
                fbOut.write(linePrevious)
                posterioris = line.strip().split()
                for i in range(0, len(posterioris), len(popList)):
                    for parental in popList:
                        classID = int(classesPerSet[setNum][parental]) - 1  # starts on 0
                        post = posterioris[classID + i]
                        fbOut.write(f'\t{post}')
                fbOut.write('\n')
            fileFB.close()
            fbOut.close()
            fbPrevious.close()
            print(f'\tDone to file {fowardBackwardWithSet}')
            print(f'\t\tRemoving the file {output}_set{previous}.fb.tsv')
            execute(f'rm {output}_set{previous}.fb.tsv')

    
    print(f'All FowardBackward was merged. Create a bim file to insert the information about chr pos and cM')
    command = f'{plink} --vcf {vcf} --cm-map {geneticMap} {chrom} --make-bed --out {output}_gen --double-id'
    execute(command)

    

    print(f'Open the {output}_gen.bim and the {output}_set{last}.fb.tsv. Creating the final FB and MSP file')
    
    
    fbFinal = open(f'{output}_AllSets.fb.tsv', 'w')
    mspFinal = open(f'{output}_AllSets.msp.tsv', 'w')
    
    fbFinal.write(f'{header}\n')
    mspFinal.write(f'{headerMSP}\n')

    fbPrevious = open(f'{output}_set{last}.fb.tsv')
    bim = open(f'{output}_gen.bim')

    markerAll = 0
    for window in SNPsPerWindowDict:
        line = fbPrevious.readline()
        for marker in range(0, SNPsPerWindowDict[window]):
            bimLine = bim.readline().split()
            if marker == 0:
                chromBim = bimLine[0]
                bpPosition= bimLine[3]
                cMPosition = float(bimLine[2])
                if cMPosition < 0:
                    cMPosition = 0
                    
        bpPositionLast = bimLine[3]
        cMPositionLast = float(bimLine[2])
        if cMPositionLast < 0:
            cMPositionLast = 0
        
        marker = marker + 1 #starts on zero
        markerAll = markerAll + marker 


        fbFinal.write(f'{chromBim}\t{bpPosition}\t{cMPosition}\t{markerAll}\t{line}')
        mspFinal.write(f'{chromBim}\t{bpPosition}\t{bpPositionLast}\t{cMPosition}\t{cMPositionLast}\t{marker}')
        
        splitLine = line.split()
        for i in range(0, len(splitLine), len(popList)):
            biggest = 0
            for j in range(1, len(popList)):
                #print(f'splitLine[{i}+{j}] ({splitLine[i+j]}) > splitLine[{i}+{biggest}] ({splitLine[i+biggest]})')
                if float(splitLine[i+j]) > float(splitLine[i+biggest]):
                    biggest = j
            #print (f'{biggest}')
            mspFinal.write(f'\t{biggest}')
            #getchar()
        mspFinal.write('\n')
            
        
        
    fbFinal.close()

    print(f'The final file was generated. Lets exclude the temp files (PLINK file and {output}_set{last}.fb.tsv)')
    command = f'rm {output}_gen.*'
    #execute(command)

    command= f'rm {output}_set{last}.fb.tsv'
    execute(command)

    fbOut.close()



def convertOutputRFMix(vcf, mapping, fowardBackward, output, setPrefix, classes, first, last, chrom, plink, geneticMap, Xlist):
    mappingDict = createMappingDict(mapping)
    classesPerSet, indPerSet, popList = createClassesPerSetDict(setPrefix, classes, mappingDict, first, last)

    print(f'Open the fowardBackward files to create a {output}.fb.tsv file')

    header = createHeader(popList)
    firstFile = True

    for setNum in range(first, last+1):
        header = addToHeader(indPerSet, popList, setNum, header)
        fowardBackwardWithSet = fowardBackward.replace('*', str(setNum))

        if firstFile:
            firstFile = False
            print(f'\tOpen the {fowardBackwardWithSet} file and create the first tempFile {output}_set{setNum}.fb.tsv')
            fileFB = open(fowardBackwardWithSet)
            fbOut = open(f'{output}_set{setNum}.fb.tsv', 'w')
            for line in fileFB:
                posterioris = line.strip().split()
                for i in range(0,len(posterioris), len(popList)):
                    for parental in popList:
                        classID = int(classesPerSet[setNum][parental]) - 1 #starts on 0
                        post = posterioris[classID+i]
                        fbOut.write(f'\t{post}')
                fbOut.write('\n')
            fileFB.close()
            fbOut.close()
            print(f'\tDone to file {fowardBackwardWithSet}')
        else:
            previous = setNum-1
            print(f'\tOpen the {fowardBackwardWithSet} file and create tempFile {output}_set{setNum}.fb.tsv')
            fileFB = open(fowardBackwardWithSet)
            fbPrevious = open(f'{output}_set{previous}.fb.tsv')
            fbOut = open(f'{output}_set{setNum}.fb.tsv', 'w')
            for line in fileFB:
                linePrevious = fbPrevious.readline().strip()
                fbOut.write(linePrevious)
                posterioris = line.strip().split()
                for i in range(0, len(posterioris), len(popList)):
                    for parental in popList:
                        classID = int(classesPerSet[setNum][parental]) - 1  # starts on 0
                        post = posterioris[classID + i]
                        fbOut.write(f'\t{post}')
                fbOut.write('\n')
            fileFB.close()
            fbOut.close()
            fbPrevious.close()
            print(f'\tDone to file {fowardBackwardWithSet}')
            print(f'\t\tRemoving the file {output}_set{previous}.fb.tsv')
            execute(f'rm {output}_set{previous}.fb.tsv')


    print(f'All FowardBackward was merged. Create a bim file to insert the information about chr pos and cM')
    command = f'{plink} --vcf {vcf} --cm-map {geneticMap} {chrom} --make-bed --out {output}_gen --double-id'
    execute(command)

    print(f'Open the {output}_gen.bim and the {output}_set{last}.fb.tsv')
    fbFinal = open(f'{output}_AllSets.fb.tsv', 'w')
    fbFinal.write(f'{header}\n')

    fbPrevious = open(f'{output}_set{last}.fb.tsv')
    bim = open(f'{output}_gen.bim')

    marker = 0
    for line in fbPrevious:
        bimLine = bim.readline().split()
        marker = marker + 1

        chromBim = bimLine[0]
        bpPosition= bimLine[3]
        cMPosition = float(bimLine[2])
        if cMPosition < 0:
            cMPosition = 0


        fbFinal.write(f'{chromBim}\t{cMPosition}\t{bpPosition}\t{marker}\t{line}')
    fbFinal.close()

    print(f'The final file was generated. Lets exclude the temp files (PLINK file and {output}_set{last}.fb.tsv)')
    command = f'rm {output}_gen.*'
    execute(command)

    command= f'rm {output}_set{last}.fb.tsv'
    execute(command)

    fbOut.close()

def execute(command):
    print(f'\t\t\tRun: {command}')
    os.system(f'{command}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='RFMix1 to RFMix2')

    required = parser.add_argument_group("Required arguments")
    required.add_argument('-v', '--vcf', help='VCF file original', required=True) #VCF

    #To PLINK
    required.add_argument('-g', '--geneticMap', help='Genetic map to generate the bim', required=True)  #Out
    required.add_argument('-c', '--chromosome', help='Number of chromosome', required=True) #VCF

    #IO
    required.add_argument('-o', '--output', help='Output prefix', required=True) #out and VCF

    #Mapping
    required.add_argument('-m', '--mapping', help='Mapping file (format: IND POP)', required=True) #out

    #Set
    required.add_argument('-s', '--set', help='Sets file with set number replaced by *', required=True)  # out and VCF
    required.add_argument('-b', '--begin', help='First set', required=True)  # out
    required.add_argument('-e', '--end', help='Last set', required=True)  # out

    #RFMix files
    required.add_argument('-C', '--classes', help='Classes file used on RFMix1 with set number replaced by *', required=True)  # out
    required.add_argument('-F', '--fowardBackward', help='ForwardBackward from RFMix1 with the set number replaced by *', required=False) #out
    #required.add_argument('-V', '--viterbi', help='Viterbi from RFMix1 with the set number replaced by *', required=False) #out
    required.add_argument('-S', '--snpPerWindow', help='SNPsPerWindow RFMix1 with the set number replaced by *', required=False, default = '') #VCF
    required.add_argument('-A', '--allelesRephased', help='allelesRephased from RFMix1 with the set number replaced by *', required=False) #VCF

    #Plink
    required = parser.add_argument_group("optional arguments")
    required.add_argument('-p', '--plink', help='Path to PLINK (default plink)', required=False, default = 'plink')
    required.add_argument('-X', '--XList', help='List of individuals created by VCF2RFMix with flag -X or --Xmen with the set number replaced by *', required=False, default = '')

    args = parser.parse_args()

    if not args.snpPerWindow:
        convertOutputRFMix(args.vcf, args.mapping, args.fowardBackward, args.output, args.set, args.classes,
                       int(args.begin), int(args.end), args.chromosome, args.plink, args.geneticMap, args.XList)
    else:
        convertOutputRFMixWindow(args.vcf, args.mapping, args.fowardBackward, args.output, args.set, args.classes,
                       int(args.begin), int(args.end), args.chromosome, args.plink, args.geneticMap, args.XList, args.snpPerWindow)
