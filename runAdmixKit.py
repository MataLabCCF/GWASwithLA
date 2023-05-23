import argparse
import gzip
import time
import os

def execute(line, logFile, run = True):
    print(f"{line}")
    if run:
        os.system(f"{line}")
    logFile.write(f"{line}\n")

def getWindow(position, dictWindow):
    window = -2
    if position < dictWindow[0]["Begin"]:
        window = 0
    elif position > dictWindow[-1]["End"]:
        window = -1
    else:
        window = -2
        for index in range(len(dictWindow)):
            # Case 1: between the begin and end of a window
            if dictWindow[index]["Begin"] <= position <= dictWindow[index]["End"]:
                window = index
                break
            # Case 2: it is between 2 windows -> give the closer window
            elif dictWindow[index]["End"] <= position <= dictWindow[index + 1]["Begin"]:
                if (dictWindow[index]["End"] - position) ** 2 < (position - dictWindow[index + 1]["Begin"]):
                    window = index
                else:
                    window = index + 1
                    break
    if window == -2:
        input(f"Error: {position}")

    return window

def createLancFile(VCF, MSP, folder, name, begin, end):
    for chrom in range(begin, end):

        MSPFile = open(MSP)

        dictAnc = {}
        dictWindow = []

        lineCount = 0
        windowCount = 0

        for line in MSPFile:
            if lineCount == 0:
                print(f"Header : {line}")
            elif lineCount == 1:
                header = line.strip().split()
            else:
                split = line.strip().split()
                dictWindow.append({"Begin" : int(split[1]), "End" : int(split[2])})

                for i in range(6, len(split), 2):
                    ind = header[i]
                    if ind not in dictAnc:
                        dictAnc[ind] = {}
                    dictAnc[ind][windowCount] = f"{split[i]}{split[i+1]}"
                windowCount = windowCount + 1
            lineCount = lineCount + 1

        gzFile = False
        if ".gz" in VCF:
            VCFFile = gzip.open(VCF)
            gzFile = True
        else:
            VCFFile = open(VCF)


        dictToLancFile = {}
        headerFlag = True
        for line in VCFFile:
            if gzFile:
                line = line.decode("utf-8")

            if headerFlag:
                if "#CHROM" in line:
                    header = line.strip().split()
                    headerFlag = False
                    lineCount = 0
            else:
                lineCount = lineCount+1
                window = getWindow(split[1], dictWindow)
                numInd = 0
                for i in range(9, len(split)):
                    numInd = numInd + 1
                    ind = header[i]

                    anc = dictAnc[ind][window]

                    if ind not in dictToLancFile:
                        dictToLancFile[ind] = {}
                        dictToLancFile[ind]["last"] = anc
                        dictToLancFile[ind]["toFile"] = ""
                    else:
                        if anc != dictToLancFile[ind]["last"]:
                            dictToLancFile[ind]["toFile"] = f"{lineCount}:{dictToLancFile[ind]['last']} "
                            dictToLancFile[ind]["last"] = anc

        lancFile = open(f"{folder}/{name}_chr{chrom}.lanc", "w")
        lancFile.write(f"{lineCount} {numInd}\n")
        for i in range(9,len(header)):
            ind = header[i]
            lancFile.write(f"{dictToLancFile[ind]['toFile']}{dictToLancFile[ind]['last']}\n")
        lancFile.close()
    return f"{folder}/{name}_chr\*.lanc"





if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run ADMIX-KIT')

    requiredGeneral = parser.add_argument_group("Required arguments")
    requiredGeneral.add_argument('-i', '--input', help='VCF file with chromosome replaced by \'*\'. ', required=True)
    requiredGeneral.add_argument('-o', '--outputName', help='Name of output name', required=True)
    requiredGeneral.add_argument('-O', '--outputFolder', help='Name of output folder', required=True)
    requiredGeneral.add_argument('-b', '--begin', help='First chromosome (default = 1)', default=1, type=int, required = True)
    requiredGeneral.add_argument('-e', '--end', help='Last chromosome  (default = 22)', default=22, type=int, required = True)
    requiredGeneral.add_argument('-l', '--localAncestry', help='Local ancestry (MSP file) file with chromosome replaced by \'*\' ', required = True)

    programs = parser.add_argument_group("Programs")
    requiredGeneral.add_argument('-p', '--plink2', help='Plink2 path', default = "plink2")
    requiredGeneral.add_argument('-a', '--admix', help='Admix-kit path', default = "admix-kit")

    args = parser.parse_args()

    lancFile = createLancFile(args.input, args.localAncestry, args.outputFolder, args.outputName, args.begin, args.end)