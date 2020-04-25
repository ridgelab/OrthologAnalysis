#! /usr/bin/env python

import os
import sys
import re
import numpy as np
import math
import argparse
def parseArgs():
	'''
	Argument parsing is done.
	Required to have an input file.
	'''
	parser = argparse.ArgumentParser(description = 'Find Identical and co-tRNA codon pairing.')
	#give input file
	parser.add_argument("-i", help = "Input Files", action= "store", dest= "input", required=False)
	#output directory
	parser.add_argument("-oc",help="Output Directory File",action="store",dest="outputdir", required=False)
	#change window
	parser.add_argument("-f",help="Ribosome Footprint",action="store",dest="footprint", type=int, default=9, required=False)
	#Have program find same aa different codon
	parser.add_argument("-c",help="Co-tRNA codon pairing",action="store_true",dest="co_trna", required=False)
	#Have program find same aa
	parser.add_argument("-comb",help="Combined co-tRNA and identical codon pairing",action="store_true",dest="comb", required=False)
	args = parser.parse_args()
	if not args.input:
		print("You must supply an input file with -i")
		sys.exit()
	return args
def outputOverallFile(orthologAverages, orthologStandardDev, currentGene, orthologIsoforms, overallFile):
	'''
	Takes four arguments:
		OrthologAverages is a dictionary of the averages of all the isoforms of one gene
		OrthologStandardDev is a dictionary of the standard deviation of all the isoforms of one gene
		currentGene is the current gene in the file
		orthologIsoforms is the list of the number of isoforms for each species
	Calculates the average
	Calculates the pooled standard Deviation
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
			n = len(n)
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
def outputToFile(allCounts, orthologAverages, orthologStandardDev, outFile):
	'''
	Takes four arguments:
		allCounts is a dictionary of how many times each codon pairs for each isoform
		orthologAverages is a dictionary of the averages of each codon pairing in a species
		orthologStandardDev is a dictionary of the standard deviations of each codon pairing in a species
		outFile is the output file
	Calculates the average and standard deviation of the codons found in a species and all of it's isoforms
	Adds the average and standard deviation to the orthologAverages and orthologStandardDev dictionaries
	returns OrthologAverages and orthologStandardDev dictionaries
	'''
	for key in allCounts:
		orthologAverages[key] = orthologAverages[key] + allCounts[key]
		codonAverage = np.average(allCounts[key])
		codonAverage = round(codonAverage, 4)
		standardDeviation = np.std(allCounts[key])
		standardDeviation = round(standardDeviation, 4)
		orthologStandardDev[key].append(standardDeviation)
		outstring = "," + str(codonAverage) + "," + str(standardDeviation)
		outFile.write(outstring)
	outFile.write("\n")
	round (standardDeviation, 4)
	return (orthologAverages, orthologStandardDev)
def checkKey(checkDictionary, key):
	'''
	Takes two arguments: a dictionary and key
	Checks to make sure the key exists in the dictionary
	'''
	if key in checkDictionary:
		return True
	else:
		return False
def getPairs(seq, allCounts):
	'''
	Takes two arguments:
		a sequence of DNA
		a dictionary to add the number of pairs to for each codon
	args.comb finds the pairing for the same amino acid for identical or different codons  
	args.co_trna finds the pairing for same amino acid different codon
	the default finds the pairing for the same codon
	args.footprint is the window for the pairing
	'''
	footprint = args.footprint
	pairs = {}
	createList = False
	sequence = []
	if args.comb or args.co_trna:
		from Bio.Seq import Seq
		from Bio.Alphabet import generic_dna
		seqaa = Seq(seq, generic_dna)
		aa = str(seqaa.translate())
		sequence = re.findall(".", aa)
		pairs = {'*':[],'X':[],'A':[],'R':[],'N':[],'D':[],'B':[],'C':[],'E':[],'Q':[],'Z':[],'G':[],'H':[],'I':[],'L':[],'K':[],'M':[],'F':[],'P':[],'S':[],'T':[],'W':[],'Y':[],'V':[]}
	if args.co_trna:
		codons = re.findall("...",seq)
	if not args.co_trna and not args.comb:

		from itertools import product
		codons = product("ACGT", repeat=3)
		for c in codons:
			c = "".join(c)
			pairs[c] = [] 
		sequence = re.findall("...", seq)
	lastFound = dict()
	for x in range(len(sequence)):
		if checkKey(pairs, sequence[x]):
			curCodon = sequence[x]
		else:
			continue
		if not curCodon in lastFound or (x - lastFound[curCodon] >= footprint):
			lastFound[curCodon] = x
			continue
		if args.co_trna:
			if codons[x] == codons[lastFound[curCodon]]:
				continue
		pairs[curCodon].append(1)
		countPairs = 0
		lastFound[curCodon] = x
	codonset = set(sequence)
	for codon in pairs:
		#if checkKey(pairs, codon):
		pairs[codon] = sum(pairs[codon])
	if allCounts == {}:
		createList = True
	for codon in pairs:
		if createList == True:
			allCounts[codon] = []
		allCounts[codon].append(pairs[codon])
	return allCounts
def readOneFile(args):
	'''
	Reads one input file that is supplied as a parameter.
	Goes through the file one line and calls functions to get codon pairing, averages, standard deviations, and output to files
	'''

	inFile = open(args.input, 'r')
	if args.outputdir:
		outDir = args.outputdir 
	else:
		outDir = 1
	if not os.path.exists(outDir): 
		os.mkdir(outDir)
	currentDir = outDir
	currentGene = ""
	currentSpecies = ""
	newSpecies = False
	firstLine = True
	orthologAverages = {}
	orthologStandardDev = {}
	overallFile = open(currentDir + "/" + "overall_file.csv", "w")
	overallFile.write("gene")
	orthologIsoforms = []
	geneCount = 0
	nameCount = 0
	isoformCount = 0
	for line in inFile:
		line = line.strip()
		line = line.split("\t")
		gene = line[0]
		if "LOC" not in gene and "CDS" not in  gene: 
			species = line[1]
			dnaSeq = line[3]
			if currentSpecies != species:
				newSpecies = True
			if currentSpecies == species and currentGene == gene:
				sequence = []
				dnaSeq = dnaSeq.replace("*", "\t")
				dnaSeq = dnaSeq.split("\t")
				for seq in dnaSeq:
					s = "".join(seq)
					sequence = sequence + re.findall("...", s)
				dnaSeq = "".join(sequence)
				allCounts = getPairs(dnaSeq, allCounts)
				isoformCount += 1.0
			else:
				if newSpecies == True and firstLine == False:
					outstring = currentGene + "," + currentSpecies + "," + str(isoformCount)
					outFile.write(outstring)
					orthologIsoforms.append(isoformCount)	
					orthologAverages, orthologStandardDev = outputToFile(allCounts, orthologAverages, orthologStandardDev, outFile)

				allCounts = {}
				isoformCount = 1.0
				sequence = []
				dnaSeq = dnaSeq.replace("*", "\t")
				dnaSeq = dnaSeq.split("\t")
				for seq in dnaSeq:
					s = "".join(seq)
					sequence = sequence + re.findall("...", s)
				dnaSeq = "".join(sequence)
				allCounts = getPairs(dnaSeq, allCounts)
				currentSpecies = species
				if currentGene != gene:
					if orthologAverages != {}:
						outputOverallFile(orthologAverages,orthologStandardDev, currentGene, orthologIsoforms, overallFile)	
						orthologIsoforms = [] 
					if geneCount == 0 or  geneCount == 500:
						nameCount += 1
						outFileName = str(nameCount) + ".csv"
						outFile = open(currentDir + "/" + outFileName, "w")
						outFileName = ""
						geneCount = 0
						outFile.write("Gene,Species,numisoforms")
						for key in allCounts:
							outFile.write(",")
							overallFile.write(",")
							outFile.write(key)
							overallFile.write(key)
							outFile.write(",std_")
							overallFile.write(",plsd_")
							outFile.write(key)
							overallFile.write(key)
						outFile.write("\n")
					for key in allCounts:
						orthologAverages[key] = []
						orthologStandardDev[key] = []
					currentGene = gene
					geneCount += 1
				firstLine = False
	orthologIsoforms.append(isoformCount)
	outstring = currentGene + "," + currentSpecies + "," + str(isoformCount)
	outFile.write(outstring)
	orthologAverages, orthologStandardDev = outputToFile(allCounts, orthologAverages, orthologStandardDev, outFile)
	outputOverallFile(orthologAverages, orthologStandardDev, currentGene, orthologIsoforms,overallFile)
	outFile.close()
	overallFile.close()
	inFile.close() 
 
'''
Main.
'''
args = parseArgs()
readOneFile(args)
