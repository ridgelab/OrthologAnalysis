#! /usr/bin/env python

import sys
import re
import os
import math
import numpy as np
from Bio.Seq import Seq

def makeCodonList(dnaSeq):
	'''
	Takes one argument: a list of dna 
	Makes a list of codon from the list of dna
	Returns codon list
	'''
	seqCodonList = []
	dnaSeq = dnaSeq.replace("*", "\t")
	dnaSeq = dnaSeq.split("\t")
	for seq in dnaSeq:
		seq = "".join(seq)
		codons = list(re.findall("...",seq))
		seqCodonList = seqCodonList + codons
	return seqCodonList 
def makeAminoAcidList(dnaSeq):
	'''
	Takes one argument a list of DNA 
	Uses BioSeq to translate DNA list to a list of AminoAcids
	Creates a dictionary
		Keys are each possible amino acids
		Values are lists. 1 appeneded to the list to count the number of each amino acid from original translated DNA sequence
	The 1's are added together and the sum is appened to a new dictionary (aminoAcidCounts)
	Returns dictionary of amino acid counts
	'''
	combinedSeq = []
	for seq in dnaSeq:
		seq = Seq(seq)
		aminoSeq = seq.translate(gap = "N") 
		combinedSeq.append(aminoSeq)
	for aminoAcid in combinedSeq:
		if aminoAcid in aminoAcidDictionary:
			aminoAcidDictionary[aminoAcid].append(1)
	for counts in aminoAcidDictionary:
		sumOfCounts = sum(aminoAcidDictionary[counts])
		aminoAcidDictionary[counts] = [0]
		aminoAcidCounts[counts].append(sumOfCounts)
		sumOfCounts = 0
	return aminoAcidCounts
def outputOverallFile(orthologAverages, orthologStandardDev, currentGene, orthologIsoforms):
	'''
	Takes four arguments:
		OrthologAverages is a dictionary of the averages of all the isoforms of one gene
		OrthologStandardDev is a dictionary of the standard deviation of all the isoforms of one gene
		currentGene is the current gene in the file
		orthologIsoforms is the list of the number of isoforms for each species
	Calculates the average of each amino acid in the gene
	Calculates the pooled standard Deviation of the gene
		n = sample size
		s = standard deviation
		k = number of samples
		Spooled = sqrt(((n1- 1)S1^2 +...+(nk - 1)Sk^2)/(n1 +...+nk-k))
	Writes the average and pooled standard deviation to a file
	'''
	overallFile.write(currentGene)
	count = 0
	k = 0
	plSd = 0
	sumOfNs = 0
	average = 0
	for key in orthologAverages:
		if key != "X" and key != "aminoAcidsAvoided":
			overallFile.write(",")
			isoformsStart = 0
			isoformsEnd = 0
			average = np.average(orthologAverages[key]) 
			average = round(average, 4)
			for length in orthologStandardDev[key]:
				n = orthologAverages[key]
				isoformsEnd = isoformsStart + int(orthologIsoforms[count]) 
				if isoformsEnd != len(orthologAverages[key]):
					n = n[isoformsStart:isoformsEnd]

				else:
					n = n[isoformsStart:]
				n = sum(n)
				sd = orthologStandardDev[key]
				sd = sd[count]
				sumOfNs = sumOfNs + n
				n = n - 1
				sd = sd * sd
				nTimesSd = n * sd
				plSd = plSd + nTimesSd 
				count += 1
				k += 1
				isoformsStart = isoformsEnd
			count = 0
			k = sumOfNs - k
			sumOfNs = 0

			if k != 0:
	
				plSd = plSd / k
				if plSd >= 0.0:	
					plSd = math.sqrt(plSd)
			if plSd == -0.0:
				plSd = 0.0
			k = 0
			plSd = round(plSd, 4)
			outputString = str(average) + "," + str(plSd)

			overallFile.write(outputString)

	
			plSd = 0
		
		else:

			if key != "X":
				average = np.average(orthologAverages["aminoAcidsAvoided"]) 
				average = round(average, 4)
				plSd = np.std(orthologAverages["aminoAcidsAvoided"])
				plSd = round(plSd, 4)
				outputString = "," + str(average) + "," + str(plSd)

				overallFile.write(outputString)

	overallFile.write("\n")
def outputToFile(aminoAcidCounts, numIsoforms, orthologStandardDev, orthologAverages):
	'''
	Takes four arguments:
		aminoAcidCounts dictionary of the counts of each amino acid in a species
		numIsoforms the number of isoform
		orthologStandardDev dictionary of standard deviations of each species 
		orthologAverages dictionary of the average of amino acids used in each species 
	Calculates the average of each amino acid and standard deviation of a species 
	Appends the averages and standard deviations to orthologAverages and orthologStandardDev
	Returns orthologStandardDev and orthologAverages 
	'''
	aminoAcidsAvoided = 0
	outstring = ""
	outFile.write(currentGene) 
	outFile.write(",")
	outFile.write(currentSpecies)
	outFile.write(",")
	numIsoforms = str(numIsoforms)
	outFile.write(numIsoforms)
	standardDeviation = 0.0
	aminoAcidAverage = 0.0
	for key in aminoAcidCounts:
		orthologAverages[key].extend(aminoAcidCounts[key])
		if aminoAcidCounts[key] != []:
			aminoAcidAverage = np.average(aminoAcidCounts[key])
			aminoAcidAverage = round(aminoAcidAverage, 4)
			standardDeviation = np.std(aminoAcidCounts[key])
			standardDeviation = round(standardDeviation, 4)
			orthologStandardDev[key].append(standardDeviation)
		else:
			aminoAcidAverage = 0
			standardDeviation = 0.0
			orthologStandardDev[key].append(standardDeviation)
		if aminoAcidAverage == 0 and key != 'X':
			aminoAcidsAvoided += 1
		if key != 'X':
			outstring = "," + str(aminoAcidAverage) + "," + str(standardDeviation)
			outFile.write(outstring)
			
	outstring = "," + str(aminoAcidsAvoided) + "\n"
	outFile.write(outstring)
	orthologAverages["aminoAcidsAvoided"].append(aminoAcidsAvoided)
	standardDeviation = np.std(aminoAcidsAvoided)
	orthologStandardDev["aminoAcidsAvoided"].append(standardDeviation)
	return (orthologStandardDev, orthologAverages)	
def outputFirstLine():
	'''
	Outputs the first line of the file
	'''	
	outFile.write("Gene,Species,numisoforms,")
	outFile.write('*,*std,A,A_std,R,R_std,N,N_std,D,D_std,B,B_std,C,C_std,E,E_std,Q,Q_std,Z,Z_std,G,G_std,H,H_std,I,I_std,L,L_std,K,K_std,M,M_std,F,F_std,P,P_std,S,S_std,T,T_std,W,W_std,Y,Y_std,V,V_std,AminoAcidsAvoided') 	
	outFile.write("\n")
	'''
	Main
	makes a directory for output files
	loops one line at a time through the file
	When there is a new species for the same gene it outputs to a file the average and standard deviation
	For each new gene it outputs the average and pooled standard deviation for the gene
	'''

inFile = open(sys.argv[1], "r")
outDir = sys.argv[2] 
if not os.path.exists(outDir): 
	os.mkdir(outDir)
currentDir = outDir
currentGene = ""
currentSpecies = ""
newSpecies = False
firstLine = True
overallFile = open(currentDir + "/" + "overall_file.csv", "w")
overallFile.write("gene,*,plstd_*,A,A_plstd,R,R_plstd,N,N_plstd,D,D_plstd,B,B_plstd,C,C_plstd,E,E_plstd,Q,Q_plstd,Z,Z_plstd,G,G_plstd,H,H_plstd,I,I_plstd,L,L_plstd,K,K_plstd,M,M_plstd,F,F_plstd,P,P_plstd,S,S_plstd,T,T_plstd,W,W_plstd,Y,Y_plstd,V,V_plstd,aminoAcidsAvoided,aminoAcidsAvoided_plsd")	
overallFile.write("\n")
orthologStandardDev = {}
geneCount = 0
nameCount = 1
outFile = open(currentDir + "/1.csv", "w")
outFile.write("Gene,Species,numisoforms,")
outFile.write('*,*std,A,A_std,R,R_std,N,N_std,D,D_std,B,B_std,C,C_std,E,E_std,Q,Q_std,Z,Z_std,G,G_std,H,H_std,I,I_std,L,L_std,K,K_std,M,M_std,F,F_std,P,P_std,S,S_std,T,T_std,W,W_std,Y,Y_std,V,V_std,AminoAcidsAvoided') 	
outFile.write("\n")
for line in inFile:
	line = line.strip()
	line = line.split("\t")
	gene = line[0]
	if "LOC" not in gene and "CDS" not in  gene: 
		species = line[1]
		dnaSeq = line[3]
		dnaSeq = makeCodonList(dnaSeq)
		if currentSpecies != species:
			newSpecies = True
		if currentSpecies == species and currentGene == gene:
			aminoAcidCounts = makeAminoAcidList(dnaSeq)
			isoformCount += 1.0
		else:
			if newSpecies == True and firstLine == False:
				orthologStandardDev, orthologAverages = outputToFile(aminoAcidCounts, isoformCount, orthologStandardDev, orthologAverages)
				isoformList.append(isoformCount) 
			isoformCount = 1
			aminoAcidDictionary = {'*':[0],'X':[0],'A':[0],'R':[0],'N':[0],'D':[0],'B':[0],'C':[0],'E':[0],'Q':[0],'Z':[0],'G':[0],'H':[0],'I':[0],'L':[0],'K':[0],'M':[0],'F':[0],'P':[0],'S':[0],'T':[0],'W':[0],'Y':[0],'V':[0]}
			aminoAcidCounts = {'*':[],'X':[],'A':[],'R':[],'N':[],'D':[],'B':[],'C':[],'E':[],'Q':[],'Z':[],'G':[],'H':[],'I':[],'L':[],'K':[],'M':[],'F':[],'P':[],'S':[],'T':[],'W':[],'Y':[],'V':[]}
			aminoAcidCounts = makeAminoAcidList(dnaSeq)
			currentSpecies = species
			if currentGene != gene:
				if orthologStandardDev != {}:
					outputOverallFile(orthologAverages, orthologStandardDev, currentGene, isoformList)
				orthologAverages = {'*':[],'X':[],'A':[],'R':[],'N':[],'D':[],'B':[],'C':[],'E':[],'Q':[],'Z':[],'G':[],'H':[],'I':[],'L':[],'K':[],'M':[],'F':[],'P':[],'S':[],'T':[],'W':[],'Y':[],'V':[],'aminoAcidsAvoided':[]}
				orthologStandardDev = {'*':[],'X':[],'A':[],'R':[],'N':[],'D':[],'B':[],'C':[],'E':[],'Q':[],'Z':[],'G':[],'H':[],'I':[],'L':[],'K':[],'M':[],'F':[],'P':[],'S':[],'T':[],'W':[],'Y':[],'V':[], 'aminoAcidsAvoided':[]}
				isoformList = []
				if geneCount == 500:
					nameCount += 1
					outFileName = str(nameCount) + ".csv"
					outFile = open(currentDir + "/" + outFileName, "w")
					outFileName = ""
					outputFirstLine()
					geneCount = 0
				currentGene = gene
				geneCount += 1
			firstLine = False

orthologStandardDev, orthologAverages = outputToFile(aminoAcidCounts, isoformCount, orthologStandardDev, orthologAverages)
isoformList.append(isoformCount) 
outputOverallFile(orthologAverages, orthologStandardDev, currentGene, isoformList)
outFile.close()
