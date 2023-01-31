#!/usr/bin/env python
# coding: utf-8

# Copyright (C) <2023>  <The Ohio State University>       
# 
# This program is free software: you can redistribute it and/or modify                              
# it under the terms of the GNU General Public License as published by 
# the Free Software Foundation, either version 3 of the License, or    
# (at your option) any later version.                                                                                       
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of           
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       
# GNU General Public License for more details.                                                                             
# 
# You should have received a copy of the GNU General Public License 
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# This script is used to go through all sequences that were used for the analysis and generate a csv file with all the relevant information to be included as a supplemental table with the manuscript.

# In[24]:


import csv


# Define the groups Ben divided the species into based on the phylogenetic tree.

# In[25]:


BensGroups=['N627','N647','N688','N825','N855','N907','N916','N927','N1144','N1278','N1311','N1325','N1394','N1432','N1441']
BensNames=['unknown','Escherichia','Mixta','Xenorhabdus','Cronobacter','Photorhabdus','Proteus','Frischella','Plesiomonas','Yersinia','Serratia','Providencia I','Budvicia','Pectobacterium','Dickeya','Providencia II']


# In[26]:


groupIndices = dict()
for i in range(len(BensGroups)):
    alignmentFile = open('classes/clades/' + BensGroups[i] + '.aln')
    for line in alignmentFile:
        if len(line) > 20:
            accessionStart = line.find("_GC")
            if accessionStart >= 0:
                accessionEnd = line.find("_", accessionStart+1)
                if accessionEnd >= 0:
                    accessionEnd = line.find("_", accessionEnd+1)
                if accessionEnd >= 0:
                    accessionEnd = line.find("_", accessionEnd+1)
                if accessionEnd >= 0:
                    accessionPlusIndex = line[accessionStart+1:accessionEnd]
                    groupIndices[accessionPlusIndex]=i


# Read the file with all sequences that were used in the analysis

# In[27]:


allSequencesFile = open('AllToBeFolded.fa', 'r')
allSequencesLines = allSequencesFile.readlines()
allSequencesFile.close()


# Go through all sequences, interpret their defline, and add to a table

# In[28]:


sequence = ''
supplementaryTable = [['Accession','Species','16S index','group index','group name','Sequence']]
for line in allSequencesLines:
    cleanLine = line.rstrip()
    if cleanLine[0] == '>':
        cleanLine = cleanLine[1:] # remove leading ">"
        if sequence != '':
            groupIndex = 0
            key = accession + "_" + str(index)
            if key in groupIndices:
                groupIndex = groupIndices[key]+1
            supplementaryTable.append([accession,species,index,groupIndex,BensNames[groupIndex],sequence])
        sequence = ''
        sequenceID = cleanLine.split()[0]
        sequenceIDs = sequenceID.split('_')
        index = sequenceIDs[-1]
        accession = sequenceIDs[-3] + '_' + sequenceIDs[-2]
        if accession[0:2]!='GC':
            print('Error in parsing ' + cleanLine + '!\n')
        species = ' '.join(sequenceIDs[:-3])
    else:
        sequence += cleanLine.replace('@','N')
# output the last one
if sequence != '':
    groupIndex = 0
    key = accession + "_" + str(index)
    if key in groupIndices:
        groupIndex = groupIndices[key]+1
    supplementaryTable.append([accession,species,index,groupIndex,BensNames[groupIndex],sequence])


# output the data into a file

# In[29]:


with open('AllUsedSequences.csv', 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile, dialect='excel')
    csvwriter.writerows(supplementaryTable)


# In[ ]:




