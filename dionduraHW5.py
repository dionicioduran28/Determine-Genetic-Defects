from itertools import islice

def readSNPedia(filename):
    file = open(filename, encoding='utf-8')
    snpData = []
    
    #iterate through snpPedia file 
    for line in islice(file,1,None):
        splitLine = line.split(",")
        for i in splitLine:
            #print (i)
            #add line to snpData if the first characters are rs or " as those
            #are the only valid lines 
            if (i[0] == '"' or i[0:2] == 'rs' ):            
                snpData.append(splitLine)

    #create new list that is clean in that all lists within list are uniform
    #begin with the snp and end with the nucleotide pair 
    cleanSnpData = []
    for i in range(len(snpData)-1):
        #if the first two letters of the list are rs and same for the next 
        #add that line to clean SNP
        
        if (snpData[i][0][0:2] == 'rs' and snpData[i+1][0][0] == 'rs'):
            cleanSnpData.append(snpData[i])
            
        #if the first two letters of the list are rs and " for the next 
        #add the concatnation of the two lines to cleanSNP
        if (snpData[i][0][0:2] == 'rs' and snpData[i+1][0][0] == '"'):
            cleanSnpData.append(snpData[i] + snpData[i+1])
    
    omitList = ['common in complete genomics' ,'common in clinvar','common on affy axiom data','common form', 'common/normal','normal','Normal risk of Basal Cell Carcinoma.', 'Normal risk of Atopic Dermatitis.', 'normal fat metabolism; common in clinvar', 'Common/Normal', 'Normal', 'can be ignored for now', 'None', 'Normal carbamazepine sensitivity; common in clinvar', 'most common in all populations', 'Normal IGF-1 level', 'normal risk', 'Normal perception of pain', "normal risk of Sjogren's syndrome", 'Normal risk of drug induced liver injury with diclofenac', 'Normal risk for autoimmune disorders', 'normal risk for infertility in Chinese men', 'normal risk of malaria', 'Normal risk for development of MGUSl', 'normal/common in clinvar', 'Normal risk for meningioma', 'Normal risk of bipolar disorder or schizophrenia', 'common', '2x risk','1.2x risk', '1.3x risk', '1.5x risk', '1.8x risk', '1.1x risk', 'reversed normal','3x risk', '> 1.48x risk', '1.18x risk', '1.6x risk' ]
    
    #create snpDict 
    snpDict = {} 
    #iterate through clean data making keys the snp and pair and value 
    #the description
    for i in cleanSnpData:
        
        #create placer variables in order to remove new lines from list and "
        #probably should have done this outside of this loop
        hold_1 = i[3]
        i[3] = hold_1.replace('"','')
        
        hold_2 = i[3]
        i[3] = hold_2.replace('\n','')
        
        hold = i[5]
        i[5] = hold.replace('\n','')
        
        #create a variable which holds the list of the the genome pair with
        #the parenthesis, done this way so that any nucelotide pair with incorrect 
        #length will be excluded later on in the conditionals 
        wholePair = i[5][-5:]

        #do not add any of summaries from omitList to the dictionary (maybe should do this later on in write summary but i honestly just thought this made more sense) and if wholePair does not begin with a parentheis 
        #you know the length is incorrect 
        if (i[3] not in omitList and wholePair[0] == '('):
      
            #create snpDict where the key is the snp smashed together with pair
            #and value is the summary 
            snpDict[i[0]+wholePair[1]+wholePair[3]] = i[3]      
    return snpDict
    
def read23(filename):
    file = open(filename, encoding='utf-8')
    meDict = {}
    
    #iterate through 23andMe file creating dictionary 
    for line in islice(file,20,None):
        splitLine = line.split()
        
        #add to the dictionary if snip starts with rs, making 
        #snp the key and value the nucleotide pair 
        if(splitLine[0][0:2] == 'rs'and splitLine[3] != '--'):
            meDict[splitLine[0]] = splitLine[3]
    return meDict

def writeSummary(filename, snpDict, meDict):
    
    #create the disease dict by finding the matches between 23andMe file and
    #snPedia 
    diseaseDict = {}
    for i in meDict:
        if ((i+meDict[i]) in snpDict):
            dictListDis = []
            #if there is a match then add to diseaseDict a key who is the snp
            #and value is a list containing the nucleotide pair and the summary 
            dictListDis.append(meDict[i])
            dictListDis.append(snpDict[i+meDict[i]])
            diseaseDict[i] = dictListDis
    
    #write results to file 
    file = open(filename, 'w', encoding='utf-8')
    for i in diseaseDict:
        file.write(str(i) + '\t' + str(diseaseDict[i][0]) + '\t' 
                   + str(diseaseDict[i][1]) + '\n')
 

def main():
    
    SNPedia = readSNPedia("snpedia.txt")
    meAnd23 = read23("genome.txt")
    results = writeSummary('genomeResults.txt', SNPedia, meAnd23)

main()
