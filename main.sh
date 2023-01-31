#!/bin/bash

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


# collect all the leader-trailer sequences for folding
/bin/rm AllToBeFolded.fa
for i in *.fasta
  do
    sed '2,3d' $i >> AllToBeFolded.fa
  done

# prepare them for locarna  
sed s/@/N/g AllToBeFolded.fa | sed s/-//g | sed "s/ /_/g" > AllToBeFolded_for_locarna.fas

# run RNAclust
RNAclust.pl --fasta AllToBeFolded_for_locarna.fas --dir AllToBeFolded_clust --cpu 8

# add labels to the tree
awk '{n=split($0,u,")");res="";for(i=1;i<n;i++) {res = res u[i] ")N" i;} print(res u[n]);}' AllToBeFolded_clust/tree > AllToBeFolded.tree

# make directories for output data
mkdir -f classes
mkdir -f classes/clades

# process each of the subtrees
cd classes
for i in N627 N647 N688 N825 N855 N907 N916 N927 N1144 N1278 N1311 N1325 N1394 N1432 N1441
  do
    ../process_subnode.sh $i
  done
  
