import sys, re, json, os
#Thanks
#First argument gives the path to the file that matches SRR* accession codes to metadata (species, sex, tissue, replicate number).
srr2metadata = open(sys.argv[1], 'r')

#Second argument is the path to the folder that contains the directories for each species
rootDir = sys.argv[2]

#Third argument is the name of the file that can be found in each SRR directory with TPM values for the transcripts of interest
quantFileName = sys.argv[3]

#First, create a dictionary to get from SRR* accession codes to metadata
srr2metadataDict = dict()
speciesList = []
for line in srr2metadata:
	if line[0] == '>': 
		#Lines starting with a '>' symbol give the name of the species, which is stored for later use
		species = line.strip()[1:]
		speciesList.append(species)
	else: 
		#Other lines match SRR* codes to metadata, which is stored in a dictionary
		lineList = line.strip().split(',')
		srr = lineList[0]
		attrList = re.split('\s+', lineList[1])
		srr2metadataDict[srr] = {'species':species, 'tissue':attrList[0], 'sex':attrList[1], 'rep':int(attrList[2])}

		
#For each species, find the path to the relevant directory and store it in a dictionary with species as keys
availableDirs = os.listdir(rootDir)
speciesDirDict = {sp:'' for sp in speciesList}

for species in speciesList:
	for d in availableDirs:
		#Search the directories in the root directory for one with a string that matches the species name - note it has to be a perfect match, ignoring cases
		if re.match(species, d, flags = re.IGNORECASE):
			speciesDirDict[species] = os.path.join(rootDir, d)
			break
			
def makeIdNameDict(id2nameFileName):
	#Function to make a dictionary with transcript IDs as keys and corresponding gene names as values
	outDict = {}
	with open(id2nameFileName, 'r') as f:
		for line in f:
			lineList = line.strip().split('\t')
			outDict[lineList[1]] = lineList[0]
	return outDict
	
def makeTranTpmDict(quantFileName):
	#Function to make a dictionary with transcript IDs as keys and corresponding TPMs as values
	outDict = {}
	with open(quantFileName, 'r') as quantFile:
		for line in quantFile:
			lineList = re.split('\s+', line.strip())
			outDict[lineList[0]] = float(lineList[2])
	return outDict
		
def makeGeneTpmDict(tranTpmDict, id2nameDict):
	#Function that takes a transcript ID/TPM dictionary and creates a dictionary with gene names as keys and total gene TPMs as values
	outDict = {}
	for k in tranTpmDict.keys():
		geneName = id2nameDict[k]
		if geneName not in outDict.keys():
			outDict[geneName] = tranTpmDict[k]
		else:
			outDict[geneName] += tranTpmDict[k]
	return outDict
	
def getTotals(speciesDict):
	#Function that totals TPM values for genes across sexes, tissues, and replicates and stores the results in a dictionary
	outDict = {}
	for sex in speciesDict.keys():
		for tissue in speciesDict[sex].keys():
			for rep in speciesDict[sex][tissue].keys():
				geneDict = speciesDict[sex][tissue][rep]['geneDict']
				for gene in geneDict.keys():
					if gene not in outDict.keys():
						outDict[gene] = geneDict[gene]
					else:
						outDict[gene] += geneDict[gene]
	return outDict
		
#Tree-like data structure: species -> sex  -> tissue -> replicate -> {{transcripts:tpm}, {genes:tpm}}
#or species -> total -> {genes:tpm}
dataDict = {}		

for species in speciesList:
	#Initialise the output dictionary for this species
	dataDict[species] = {}
	
	#Get a list of the contents of this species' directory and separate it into files and directories
	dirContents = os.listdir(speciesDirDict[species])
	onlyFiles = [f for f in dirContents if os.path.isfile(os.path.join(speciesDirDict[species], f))]
	onlyDirs = [f for f in dirContents if not os.path.isfile(os.path.join(speciesDirDict[species],f))]

	#Some error-catching: this relies on there being exactly one file in the species directory, which is the one we want. Any other number will terminate the program. Not very flexible, but quicker to write
	if len(onlyFiles) == 0:
		print("No files found in species directory. Cannot build gene name/transcript id dictionary. Terminating")
		sys.exit()
	elif len(onlyFiles) == 1:
		id2nameDict = makeIdNameDict(os.path.join(speciesDirDict[species], onlyFiles[0]))
	else:
		print("Too many files found in species directory. Please remove unneeded files and try again. Terminating")
		sys.exit()
	
	#We now go into each subdirectory and get the quantification file we want
	for d in onlyDirs:
		#Check the directory name has the kind of accession code we want and pull it out
		try:
			srrSpan = list(re.search('SRR\d+', d).span())
		except AttributeError:
			print("The directory %s does not appear to be the kind we want. Skipping to next directory")
			continue
		srrCode = d[srrSpan[0]:srrSpan[1]]
		
		#Use the SRR code to get the relevant metadata from the dictionary we created earlier
		metadata = srr2metadataDict[srrCode]
		
		#Get the full path for the quantification file we are interested in (NB: relies on all of the directories having identically named files)
		quantFilePath = os.path.join(speciesDirDict[species], d, quantFileName)
		
		#Now build up the tree-like dictionary structure
		if metadata['sex'] not in dataDict[species].keys():
			dataDict[species][metadata['sex']] = {}
			
		if metadata['tissue'] not in dataDict[species][metadata['sex']].keys():
			dataDict[species][metadata['sex']][metadata['tissue']] = {}
			
		if metadata['rep'] not in dataDict[species][metadata['sex']][metadata['tissue']].keys():
			dataDict[species][metadata['sex']][metadata['tissue']][metadata['rep']] = {'geneDict':{}}
			
		#Use the previously-defined functions to get the transcript- and gene-level TPM values
		tranDict = makeTranTpmDict(quantFilePath)	
		dataDict[species][metadata['sex']][metadata['tissue']][metadata['rep']]['geneDict'] = makeGeneTpmDict(tranDict, id2nameDict)
		
	#Get the totals across sexes, tissues, and reps for this species	
	dataDict[species]['total'] = {'geneDict':getTotals(dataDict[species])}
		
#Pretty print the whole lot in JSON format
print(json.dumps(dataDict, sort_keys = True, indent = 4, separators = (',', ': ')))


