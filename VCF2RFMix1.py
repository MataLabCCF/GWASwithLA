# -*- coding: utf-8 -*-
"""
Created on Tue Jul 27 10:41:29 2021

@author: PEIXOTT
"""
import os
import sys
import vcfpy
import argparse

def getchar():
    c = sys.stdin.read(1)

def printHaploid(ind, indClass, firstInd):
    lineReturn = ''
    if not firstInd:
        lineReturn= ' '
    if ind not in indClass:
        lineReturn = lineReturn+"0"
    else:
        lineReturn = lineReturn+f'{indClass[ind]}'

    return lineReturn, False


def printDiploid(ind, indClass, firstInd):
    lineReturn = ''
    if not firstInd:
        lineReturn = ' '
    if ind not in indClass:
        print(f"{ind} has no class ")
        lineReturn = lineReturn + "0 0"
    else:
        print(f"{ind} is class {indClass[ind]}")
        lineReturn = lineReturn + f'{indClass[ind]} {indClass[ind]}'

    return lineReturn, False

def printClassFile(fileName, indClass, chrX, callList, chrXFileName):

    firstInd = True
    classes = open(f'{fileName}', 'w')
    if chrX:
        XList = open(f'{chrXFileName}', 'w')
    for ind in callList:
        if chrX:
            GT = callList[ind].data.get('GT').split('|')
            if len(GT) == 1:
                line, firstInd = printHaploid(ind, indClass, firstInd)
                XList.write(f'{ind}\tH\n')
            else:
                line, firstInd = printDiploid(ind, indClass, firstInd)
                XList.write(f'{ind}\tD\n')
        else:
            line, firstInd = printDiploid(ind, indClass, firstInd)
        classes.write(line)
        
    classes.close()

def createBin(vcf, genmap, output, chrom, plink):
    if genmap[-3:] == '.gz':
        command = f'cp {genmap} ./GenMap{chrom}.gz | gunzip ./GenMap{chrom}.gz'
        execute(command)
        genmap = f'./GenMap{chrom}'
    command = f'{plink} --vcf {vcf} --cm-map {genmap} {chrom} --make-bed --out {output}_gen --double-id'
    execute(command)
    return f'{output}_gen.bim'

def execute(command):
    print(command + "\n\n")
    os.system(command)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='VCF to RFMix1')

    required = parser.add_argument_group("Required arguments")
    # required.add_argument('-s', '--steps', help='File with the order to perform the analisys', required=True)
    required.add_argument('-v', '--vcf', help='VCF file to be converted',
                          required=True)
    required.add_argument('-c', '--correspondence', help='Correspondence between parental and IDs', required=True)
    required.add_argument('-C', '--chromosome', help='It is self-explanatory ', required=True)
    required.add_argument('-o', '--output', help='Output name file', required=True)
    required.add_argument('-p', '--plink', help='PLINK path', required=True)
    required.add_argument('-m', '--map', help='Genetic map', required=False)

    optional = parser.add_argument_group("Optional arguments")
    required.add_argument('-X', '--Xmen', help='This flag allows the output of haploid individuals', required=False, default=False, action='store_true')

    args = parser.parse_args()

    chrX = args.Xmen
    
    indClass = {}
    mapDict = {}
    popDict = {}
    numParental = 0
    
    correspondence = open(f'{args.correspondence}')
    
    for line in correspondence:
        split = line.strip().split()
        if split[1] not in popDict:
            numParental = numParental+1
            popDict[split[1]] = numParental
        indClass[split[0]] = popDict[split[1]]
    
    bim = createBin(args.vcf, args.map, args.output, args.chromosome, args.plink)
    
    print(f'Abrindo o {bim}')
    mapFile = open(f'{bim}')
    header = False
    for line in mapFile:
        if header:
            header = False
        else:
            split = line.strip().split()
            if(float(split[2]) < 0):
                split[2] = 0.0
            mapDict[int(split[3])] = float(split[2])
    
    correspondence.close()
    mapFile.close()
    
    count = 0
    reader = vcfpy.Reader.from_path(f'{args.vcf}')
    alleles = open(f'{args.output}_alleles', 'w')
    location = open(f'{args.output}_location', 'w')


    classPrint = False
    for record in reader:
        location.write(f'{mapDict[record.POS]}\n')
        callList = record.call_for_sample

        if not classPrint:
            printClassFile (f'{args.output}_classes', indClass, chrX, callList, f'{args.output}_X_Inds.txt')
            classPrint = True

        for ind in callList:
            GT = callList[ind].data.get('GT').split('|')
            if chrX:
                if len(GT) == 1:
                    alleles.write(f'{GT[0]}')
                else:
                    alleles.write(f'{GT[0]}{GT[1]}')
            else:
                alleles.write(f'{GT[0]}{GT[1]}')
        alleles.write('\n')
    print(f'{count} unmapped (0 cM)')
    alleles.close()
    location.close()
