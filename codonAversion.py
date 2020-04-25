#! /usr/bin/env python
import sys
import os
import re
import numpy as np
import math

def checkKey(checkDictionary, key):
	'''
	Takes two arguments: a dictionary and key 
	Checks to make sure the key exists in the dictionary
	'''
	if key in checkDictionary:
		return True
	else:
		return False

def matchCodonCount(seqCodonList, allCounts):
	'''
	Takes two argument: a list of codons and a dictionary of the counts of each codon
	Creates a dictionary of all possible codons
	Counts the number of each codon in the list of codons given 
	Adds the count to allCounts dictionary 
	Returns the dictionary with the count of each codon
	'''
	total = 0
	codonsComb = {}
	createList = False
	from itertools import product
	codons = product("ACGT", repeat=3)
	for c in codons:
		c = "".join(c)
		codonsComb[c] = 0

	codonError = ""
	for codon in seqCodonList:	
		if checkKey(codonsComb, codon):
			total = codonsComb[codon] + 1
			codonsComb[codon] = total

		else:
			codonError = codonError + codon
	if allCounts == {}:
		createList = True
	for codon in codonsComb:
		if createList == True:
			allCounts[codon] = []
		allCounts[codon].append(codonsComb[codon])
	return allCounts
def outputOverallFileFirstLine(allCounts):
	'''
	Takes one argument: a dictionary with all possible codons as keys
	Outputs the first line to the overall file
	'''
	for codon in allCounts:
		overallFile.write(codon)
		overallFile.write(",plstd_")
		overallFile.write(codon)
		overallFile.write(",")
	overallFile.write("Codons Avoided")
	overallFile.write (",CodonsAvoided_plstd")
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
	overallFile.write("\n")
	overallFile.write(currentGene)
	count = 0
	k = 0
	plSd = 0
	sumOfAverages = 0
	average = 0
	for key in orthologAverages:
	#if len(orthologAverages[key]) > 9:
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
			sumOfAverages = sumOfAverages + n
			n = n - 1
			sd = sd * sd
			nTimesSd = n * sd
			plSd = plSd + nTimesSd 
			count += 1
			k += 1
			isoformsStart = isoformsEnd
		count = 0
		k = sumOfAverages - k
		sumOfAverages = 0
		if k != 0:
		
			plSd = plSd / k
			if plSd >= 0.0:	
				plSd = math.sqrt(plSd)
		if plSd == -0.0 or plSd <= 0:
			plSd = 0.0
		k = 0
		plSd = round(plSd, 4)
		outputString = str(average) + "," + str(plSd)
		overallFile.write(outputString)
		plSd = 0
	#	else:
	#		average = np.average(orthologAverages["codonsAvoided"]) 
	#		average = round(average, 4)
	#		plSd = np.std(orthologAverages["codonsAvoided"])
	#		plSd = round(plSd, 4)
	#		outputString = "," + str(average) + "," + str(plSd)
	#		overallFile.write(outputString)
def makeCodonList(dnaSeq):
	'''
	Takes one argument a dnaSeq
	Splits the dnaSeq into codons 
	Returns the list of codons
	'''
	seqCodonList = []
	dnaSeq = dnaSeq.replace("*", "\t")
	dnaSeq = dnaSeq.split("\t")
	for seq in dnaSeq:
		seq = "".join(seq)
		codons = list(re.findall("...",seq))
		seqCodonList = seqCodonList + codons
	return seqCodonList 
def outputToFile(allCounts, numIsoforms, orthologAverages, orthologStandardDev):
	'''
	Takes four arguments:
		orthologAverages is a dictionary of the averages of each codon in a species
		orthologStandardDev is a dictionary of the standard deviations of each codon in a species
		currentGene is the gene being calculated 
		orthologIsoforms is a list of the number of isoforms in the species
	Calculates the average and standard deviation of the codons found in a species and all of it's isoforms
	Adds the average and standard deviation to the orthologAverages and orthologStandardDev dictionaries
	returns OrthologAverages and orthologStandardDev dictionaries
	'''
	codonsAvoided = 0
	outstring = currentGene + "," + currentSpecies + "," + str(numIsoforms)
	outFile.write(outstring)
	for key in allCounts:
		orthologAverages[key] = orthologAverages[key] + allCounts[key]
		codonAverage = np.average(allCounts[key])
		codonAverage = round(codonAverage, 4)
		standardDeviation = np.std(allCounts[key])
		standardDeviation = round(standardDeviation, 4)
		orthologStandardDev[key].append(standardDeviation)
		if codonAverage == 0:
			codonsAvoided += 1
		outstring = "," + str(codonAverage) + "," + str(standardDeviation)
		outFile.write(outstring)
	outstring = "," +str(codonsAvoided) + "\n"
	outFile.write(outstring)
	orthologAverages["codonsAvoided"].append(codonsAvoided)
	standardDeviation = np.std(codonsAvoided)
	round (standardDeviation, 4)
	orthologStandardDev["codonsAvoided"].append(standardDeviation) 
	return (orthologAverages, orthologStandardDev)
	'''
	Main
	Opens the file
	Creates a directory for output files
	Loops through one line at a time of file
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
orthologAverages = {}
orthologStandardDev = {}
writeToFile = True
overallFile = open(currentDir + "/" + "overall_file.csv", "w")
overallFile.write("gene,")
orthologIsoforms = []
geneCount = 0
nameCount = 0
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
			allCounts = matchCodonCount(dnaSeq, allCounts)
			isoformCount += 1.0
		else:
			if newSpecies == True and firstLine == False:
				orthologIsoforms.append(isoformCount)	
				orthologAverages, orthologStandardDev = outputToFile(allCounts, isoformCount, orthologAverages, orthologStandardDev)

			allCounts = {}
			isoformCount = 1.0
			allCounts = matchCodonCount(dnaSeq, allCounts)
			if writeToFile == True:
				outputOverallFileFirstLine(allCounts)	
			writeToFile = False
			currentSpecies = species
			if currentGene != gene:
				if orthologAverages != {}:
					#if len(orthologAverages["AAA"]) > 9:
					#	print(orthologAverages)
					outputOverallFile(orthologAverages,orthologStandardDev, currentGene, orthologIsoforms)	
					orthologIsoforms = [] 
				#if "/" in gene:
				#	for letter in gene:
				#		if letter == "/":
				#			gene = gene.replace(letter, "-")
				if geneCount == 0 or  geneCount == 500:
					nameCount += 1
					outFileName = str(nameCount) + ".csv"
					outFile = open(currentDir + "/" + outFileName, "w")
					outFileName = ""
					geneCount = 0
					outFile.write("Gene,Species,numisoforms,")
					for key in allCounts:
						outFile.write(key)
						outFile.write(",std_")
						outFile.write(key)
						outFile.write(",")
					outFile.write("CodonsAvoided")
					outFile.write("\n")
				for key in allCounts:
					orthologAverages[key] = []
					orthologStandardDev[key] = []
				orthologAverages["codonsAvoided"] = []
				orthologStandardDev["codonsAvoided"] = []
				currentGene = gene
				geneCount += 1
			firstLine = False		
orthologIsoforms.append(isoformCount)
orthologAverages, orthologStandardDev = outputToFile(allCounts, isoformCount, orthologAverages, orthologStandardDev)
outputOverallFile(orthologAverages, orthologStandardDev, currentGene, orthologIsoforms)
outFile.close()
