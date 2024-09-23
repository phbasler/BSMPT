#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
# Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

import pandas as pd
import argparse

'''
The parameters Type, Lambda1 to Lambda8, tanbeta, m12squared and vs should
have the label of the corresponding parameter.
'''

Type = 'yuktype'
Lambda1 = 'L1'
Lambda2 = 'L2'
Lambda3 = 'L3'
Lambda4 = 'L4'
Lambda5 = 'L5'
Lambda6 = 'L6'
Lambda7 = 'L7'
Lambda8 = 'L8'
tanbeta = 'tbeta'
m12squared = 'm12sq'
vs = 'vs'


def convert(IndexCol, InputFILE, OutputFILE, Seperator):
    df = pd.DataFrame()
    if IndexCol == 'False':
        df = pd.read_table(InputFILE, index_col=False, sep=Seperator)
    else:
        df = pd.read_table(InputFILE, index_col=int(IndexCol), sep=Seperator)

    frontcol = [
            Type, Lambda1, Lambda2, Lambda3, Lambda4, Lambda5, Lambda6,
            Lambda7, Lambda8, vs, tanbeta, m12squared
    ]

    Col = [c for c in frontcol if c in df] + [c for c in df if c not in frontcol]
    df = df[Col]

    df.to_csv(OutputFILE, index=False, sep='\t')


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--indexcol', help='Column which stores the index of your data', default='False')
parser.add_argument('-in', '--input', help='Input file')
parser.add_argument('-out', '--output', help='Output file')
parser.add_argument('-sep', '--seperator', help='Column separator of input file', default='\t')

if __name__ == "__main__":
    args = parser.parse_args()
    convert(args.indexcol, args.input, args.output, args.seperator)
