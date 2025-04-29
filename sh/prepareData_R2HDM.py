#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
# Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

import pandas as pd
import argparse

"""
The parameters Type, Lambda1 to Lambda5, tanbeta and m12squared should
have the label of the corresponding parameter.
"""

type_names = ["yuktype", "type", "Type"]
Lambda1_names = ["L1"]
Lambda2_names = ["L2"]
Lambda3_names = ["L3"]
Lambda4_names = ["L4"]
Lambda5_names = ["L5"]
tanbeta_names = ["tbeta"]
m12squared_names = ["m12sq", "m12sqr"]

def get_column_names(name_list, df):
    name = next((el for el in name_list if el in df), None)
    if name is None:
        raise ValueError("No column matching " + str(name_list) + " found in file. Check the input or add column name to the list.")
    return name

def convert(IndexCol, InputFILE, OutputFILE, Seperator):
    df = pd.DataFrame()
    if IndexCol == "False":
        df = pd.read_table(InputFILE, index_col=False, sep=Seperator)
    else:
        df = pd.read_table(InputFILE, index_col=int(IndexCol), sep=Seperator)

    Type = get_column_names(type_names, df)
    Lambda1 = get_column_names(Lambda1_names, df)
    Lambda2 = get_column_names(Lambda2_names, df)
    Lambda3 = get_column_names(Lambda3_names, df)
    Lambda4 = get_column_names(Lambda4_names, df)
    Lambda5 = get_column_names(Lambda5_names, df)
    m12squared = get_column_names(m12squared_names, df)
    tanbeta = get_column_names(tanbeta_names, df)

    frontcol = [Type, Lambda1, Lambda2, Lambda3, Lambda4, Lambda5, m12squared, tanbeta]

    Col = [c for c in frontcol if c in df] + [c for c in df if c not in frontcol]
    df = df[Col]

    df.to_csv(OutputFILE, index=False, sep="\t")


parser = argparse.ArgumentParser()
parser.add_argument(
    "-i",
    "--indexcol",
    help="Column which stores the index of your data",
    default="False",
)
parser.add_argument("-in", "--input", help="Input file")
parser.add_argument("-out", "--output", help="Output file")
parser.add_argument(
    "-sep", "--seperator", help="Column separator of input file", default="\t"
)

if __name__ == "__main__":
    args = parser.parse_args()
    convert(args.indexcol, args.input, args.output, args.seperator)
