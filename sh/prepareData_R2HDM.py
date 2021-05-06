# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

import pandas as pd
import numpy as np

####### The parameters Type,Lambda1 to Lambda5, tanbeta and m12squared should have the label of the corresponding 
####### parameter. If your first column is an index column you should set HasIndexCol=True
####### otherwise set HasIndexCol=False . With Seperator you have to tell which seperator your data file
#######  is using (e.g. , \t or space). Your InputFILE will then be saved to OutputFILE


HasIndexCol=True
Seperator='\t'
InputFILE='../example/R2HDM_Input.dat'
OutputFILE='R2HDM_Ordered.dat'
Type='yuktype'
Lambda1='L1'
Lambda2='L2'
Lambda3='L3'
Lambda4='L4'
Lambda5='L5'
tanbeta='tbeta'
m12squared='m12sq'


df=pd.DataFrame()
if HasIndexCol:
	df=pd.read_table(InputFILE,sep=Seperator)
else:
	df=pd.read_table(InputFILE,index_col=False,sep=Seperator)

frontcol=[Type,Lambda1,Lambda2,Lambda3,Lambda4,Lambda5,m12squared,tanbeta]

Col = [c for c in frontcol if c in df] + [c for c in df if c not in frontcol]
df=df[Col]

df.to_csv(OutputFILE,index=False,sep='\t')
