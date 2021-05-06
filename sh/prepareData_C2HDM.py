# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

import pandas as pd
import numpy as np


####### The parameters Type,Lambda1 to Lambda4,re_Lambda5,im_Lambda5, tanbeta and re_m12squared should have the label
####### of the corresponding parameter. If your first column is an index column you should set HasIndexCol=True
####### otherwise set HasIndexCol=False . With Seperator you have to tell which seperator your data file 
####### is using (e.g. , \t or space). Your InputFILE will then be saved to OutputFILE.

HasIndexCol=False
Seperator='\t'
InputFILE='../example/C2HDM_Input.dat'
OutputFILE='C2HDM_Ordered.dat'
Type='Type'
Lambda1='L1'
Lambda2='L2'
Lambda3='L3'
Lambda4='L4'
re_Lambda5='re_L5'
im_Lambda5='im_L5'
tanbeta='tbeta'
re_m12squared='re_m12sq'


df=pd.read_table(InputFILE,index_col=HasIndexCol,sep=Seperator)

frontcol=[Type,Lambda1,Lambda2,Lambda3,Lambda4,re_Lambda5,im_Lambda5,re_m12squared,tanbeta]

Col = [c for c in frontcol if c in df] + [c for c in df if c not in frontcol]
df=df[Col]

df.to_csv(OutputFILE,index=False,sep='\t')
