#!/usr/bin/env python

######
## python table_prep.py otu_table.txt
###### 

import pandas as pd
import sys

ttt = pd.read_csv(sys.argv[1],skiprows=1,header=None,index_col=False,sep="\t",nrows=25)

ttt.iat[0,0] = ""
ttt_OTU = ttt.drop(ttt.columns[len(ttt.columns)-1], axis=1)

ttt_Tax1 = ttt.drop(ttt.columns[1:len(ttt.columns)-1], axis=1) # keep only the last column with taxonomy data
ttt_Tax2 = ttt_Tax1.drop(ttt_Tax1.index[0]) # delete first row that now just says 'taxonomy'
ttt_Tax3= ttt_Tax2[ttt_Tax2.columns[len(ttt_Tax2.columns)-1]].str.split(';', n=7, expand=True) # split taxonomic levels 
ttt_Tax3.columns = ['Kingdom','Phylum','Class','Order','Family','Genus','Species'] # Name taxonomic levels
ttt_Tax4 = ttt_Tax3.replace({'k__':'',' p__':'',' c__':'',' o__':'',' f__':'',' g__':'',' s__':''}, regex=True)

ttt.to_csv('ttt.csv',sep=",",header=False) # print sample data to add to final MiSeq sample sheet
ttt_OTU.to_csv('table_otu.csv',sep=",",header=False,index=False) # print sample data to add to final MiSeq sample sheet
ttt_Tax4.to_csv('table_tax.csv',sep=",",header=True,index=False) # print sample data to add to final MiSeq sample sheet