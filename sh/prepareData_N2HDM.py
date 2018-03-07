import pandas as pd
import numpy as np

####### The parameters Type,Lambda1 to Lambda8, tanbeta, m12squared and v_s should have the label of the 
####### corresponding parameter. If your first column is an index column you should set HasIndexCol=True
####### otherwise set HasIndexCol=False . With Seperator you have to tell which seperator your data file
####### is using (e.g. , \t or space). Your InputFILE will then be saved to OutputFILE.

HasIndexCol=False
Seperator='\t'
InputFILE='../example/N2HDM_Input.dat'
OutputFILE='N2HDM_Ordered.dat'
Type='p_THDMtype'
Lambda1='p_L1'
Lambda2='p_L2'
Lambda3='p_L3'
Lambda4='p_L4'
Lambda5='p_L5'
Lambda6='p_L6'
Lambda7='p_L7'
Lambda8='p_L8'
tanbeta='p_tbeta'
m12squared='p_m12sq'
vs='p_vs'


df=pd.read_table(InputFILE,index_col=HasIndexCol,sep=Seperator)

frontcol=[Type,Lambda1,Lambda2,Lambda3,Lambda4,Lambda5,Lambda6,Lambda7,Lambda8,vs,tanbeta,m12squared]

Col = [c for c in frontcol if c in df] + [c for c in df if c not in frontcol]
df=df[Col]
print(df)

df.to_csv(OutputFILE,index=False,sep='\t')
