#!/bin/bash
#
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

#
# check command line arguments
#
if [[ $# -ne 1 ]]; then
    echo "Illegal number of arguments"
    exit 1
fi
if ! [[ $1 =~ ^N[0-9]*$ ]]; then
    echo "Argument must be an N followed by a number"
    exit 2
fi
#
# remove previous data
#
/bin/rm -rf $1
#
# make new directory
#
mkdir $1
cd $1
#
# get the ids of all the organisms in the clade
#
nw_clade ../AllToBeFolded.tree $1 | grep -o "[(,]id[[:digit:]]*" | grep -o "[[:digit:]]*" | sort -n > all_org_ids.txt
#
# get the actual sequences
#
awk ' \
BEGIN { \
  n=0; \
  f=0; \
  i=1; \
  c=0; \
} \
ARGIND==1 { \
  n++; orgid[n]=$1; \
} \
ARGIND==2 && $0~/>/ { \
  c++; \
  if (orgid[i]==c) { \
    f=1; \
    i++; \
  } else { \
    f=0; \
  } \
} \
ARGIND==2 && f==1 { \
 print($0); \
}' all_org_ids.txt ../../AllToBeFolded.fa | sed s/@/N/g | sed "s/ /_/g" > all_sequences.fa
#
# run actual locarna
#
mlocarna -tgtdir $1.out --cpu `cat /proc/cpuinfo | grep processor | wc -l` all_sequences.fa
#
# copy results
#
cp -f $1.out/results/result.aln ../clades/$1.aln
cp -f $1.out/results/alirna.ps ../clades/${1}_alirna.ps
cp -f $1.out/results/aln.ps ../clades/${1}_aln.ps
