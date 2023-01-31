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
# 

# Load an assembly and its annotation from NCBI, identify the 16S rRNA, and extract its leader and trailer

# In[1]:


import os
import urllib
import gzip
import regex
import pandas as pd
from Bio import SeqIO


# The following two functions are from https://dmnfarrell.github.io/bioinformatics/assemblies-genbank-python and then modified to fit our needs

# In[2]:


def get_assembly_summary(id):
    """Get esummary for an entrez id"""
    from Bio import Entrez
    esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
    esummary_record = Entrez.read(esummary_handle)
    return esummary_record


# In[3]:


def get_assemblies(term, download=True):
    """Download genbank assemblies for a given search term.
    Args:
        term: search term, usually organism name
        download: whether to download the results
    """

    from Bio import Entrez
    Entrez.email = "bundschuh.2@osu.edu"
    handle = Entrez.esearch(db="assembly", term=term, retmax='200')
    record = Entrez.read(handle)
    ids = record['IdList']
    filenames = []
    for id in ids:
        #get summary
        summary = get_assembly_summary(id)
        #get ftp link
        url = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
        if url == '':
            continue
        label = os.path.basename(url)
        #get the fasta link
        link = os.path.join(url,label+'_genomic.fna.gz')
        if download == True:
            #download link
            urllib.request.urlretrieve(link, f'{label}.fna.gz')
        #get the feature table link
        link = os.path.join(url,label+'_feature_table.txt.gz')
        if download == True:
            #download link
            urllib.request.urlretrieve(link, f'{label}_ft.txt.gz')
        # remember file names
        filenames.append([summary['DocumentSummarySet']['DocumentSummary'][0]['SpeciesName'],summary['DocumentSummarySet']['DocumentSummary'][0]['SpeciesTaxid'],label+'.fna.gz',label+'_ft.txt.gz'])
    return filenames


# In[4]:


get_assemblies('GCF_003956045.1')


# In[5]:


def get16SLeaderAndTrailer(fastafile, featurefile):
    """Extract 16S rRNA leader and trailer for a downloaded organism.
    Args:
        fastafile: genome sequence file in fasta format
        featurefile: feature table file
    """
    # read feature table
    ft = pd.read_csv(featurefile, sep='\t')
    # get first 16S ribosomal RNA
    rRNA16S = ft[ft['name']=='16S ribosomal RNA']
    results = []
    for index,rRNA16S in rRNA16S.iterrows():
        start = rRNA16S['start']
        end = rRNA16S['end']
        leader=''
        trailer=''
        with gzip.open(fastafile, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if record.id == rRNA16S['genomic_accession']:
                    if rRNA16S['strand']=='+':
                        leader=record.seq[start-280:start+17]
                        preliminarytrailer=record.seq[max(0,end-100):end+200]
                    else:
                        leader=record.seq[end-18:end+279].reverse_complement()
                        preliminarytrailer=record.seq[max(0,start-200):start+100].reverse_complement()
        threeprimemotifpositions = regex.findall('(AAGTCGTAACAAGGTA){e<=2}', str(preliminarytrailer))
        if (threeprimemotifpositions == []):
            print('3\' motif not found for ' + fastafile)
        elif (len(threeprimemotifpositions)>1):
            print('multiple 3\' motifs found for ' + fastafile)
        else:
            threeprimemotifposition = regex.search('(?b)(AAGTCGTAACAAGGTA){e<=2}', str(preliminarytrailer)).end()
            trailer = preliminarytrailer[threeprimemotifposition+35:threeprimemotifposition+92]
            if len(leader)+len(trailer) > 350:
                results.append(str(leader) + "@@@@@@" + str(trailer))
    return results


# In[6]:


get16SLeaderAndTrailer('GCF_003956045.1_ASM395604v1.fna.gz','GCF_003956045.1_ASM395604v1_ft.txt.gz')


# In[7]:


def get16SSequence(term):
    """Download genome and annotation and extract 16S rRNA leader and trailer.
    Args:
        term: search term describing organism
    """
    results = []
    # download all the genomes
    genomeinfo = get_assemblies(term)
    # process all the genomes
    for genome in genomeinfo:
        rRNA16Sdata = get16SLeaderAndTrailer(genome[2],genome[3])
        if len(rRNA16Sdata) > 0:
            for i, seq in enumerate(rRNA16Sdata, start=1):
                if (seq.count('N')>5):
                    print('Too many Ns in ' + genome[0] + ' - ' + seq + '.')
                    continue
                boxCpositions = regex.findall('(TCTGTGTGGG){s<=2}',seq[:190])
                if (boxCpositions == []):
                    print('Box C not found for ' + genome[0] + ' (' + term + ') - ' + seq)
                    continue
                if (len(boxCpositions)>1):
                    print(str(len(boxCpositions)) + ' boxes C (' + ','.join(boxCpositions) + ') found for '+ genome[0] + ' (' + term + ') - using the 3\' most one.')
                boxCposition=regex.search('(?r)(TCTGTGTGGG){s<=2}',seq[:190]).start()
                if (boxCposition < 23):
                    print('Box C too far upstream for ' + genome[0] + ' (position ' + str(boxCposition) + ','+ term + ')')                    
                    continue
                boxAseq = 'N'
                boxApos = -1
                boxAseqs = regex.findall('(TGCTCTTTA){s<=2}',seq[:boxCposition])
                if (len(boxAseqs)==1):
                    boxAseq = boxAseqs[0]
                    boxApos = boxCposition - regex.search('(TGCTCTTTA){s<=2}',seq[:boxCposition]).end()
                results.append([genome[0], genome[1], term, i, boxAseq, seq[:boxCposition-23], seq[boxCposition-23:boxCposition], seq[boxCposition:boxCposition+10], seq[boxCposition+10:], boxApos])
                name = genome[0].replace(' ','_') + '_' + term + '_' + str(i)
                f = open(name + '.fasta', 'w')
                f.write('>' + name + ' (' + genome[1] + ')\n')
                f.write(boxAseq + '\n' + seq[:boxCposition-23] + '\n' + seq[boxCposition-23:boxCposition] + '\n' + seq[boxCposition:boxCposition+10] + '\n' + seq[boxCposition+10:] + '\n')
        else:
            print('No 16S information found for ' + genome[0])
        os.remove(genome[2])
        os.remove(genome[3])
    return results


# In[8]:


print(get16SSequence('GCF_001975225.1'))
print(get16SSequence('GCA_900447185.1'))
print(get16SSequence('GCF_900475855.1'))
print(get16SSequence('GCF_000400505.1'))
print(get16SSequence('GCA_006974205.1'))
print(get16SSequence('GCF_900475855.1'))
print(get16SSequence('GCF_000743055.1'))


# Now, we are going to call this function for all GTDB species representatives from the Enterobacteriaceae family of the GTDB taxonomy

# In[9]:


#all_organisms=pd.read_csv('AllGTDBEnterobacteriaceae.csv')
all_organisms=pd.read_csv('AllGTDBEnterobacteriaceaeWithNCBIType.csv')
boxAdistanceCounts = [0 for i in range(400)]
for index, row in all_organisms.iterrows():
    print(row['ID'])
    for rRNAdata in get16SSequence(row['ID']):
        boxAdistanceCounts[rRNAdata[9]+1] = boxAdistanceCounts[rRNAdata[9]+1] + 1
print(boxAdistanceCounts)


# In[ ]:




