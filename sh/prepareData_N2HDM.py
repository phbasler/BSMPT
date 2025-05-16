#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
# Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

import pandas as pd
import argparse

"""
The column labels used in the input file need to be contained in the parameter name lists.
"""

type_names = ["yuktype", "type", "Type", "p_THDMtype"]
Lambda1_names = ["L1", "p_L1"]
Lambda2_names = ["L2", "p_L2"]
Lambda3_names = ["L3", "p_L3"]
Lambda4_names = ["L4", "p_L4"]
Lambda5_names = ["L5", "p_L5"]
Lambda6_names = ["L6", "p_L6"]
Lambda7_names = ["L7", "p_L7"]
Lambda8_names = ["L8", "p_L8"]
tanbeta_names = ["tbeta", "p_tbeta"]
m12squared_names = ["m12sq", "m12sqr", "p_m12sq"]
vs_names = ["vs", "p_vs"]

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
    Lambda6 = get_column_names(Lambda6_names, df)
    Lambda7 = get_column_names(Lambda7_names, df)
    Lambda8 = get_column_names(Lambda8_names, df)
    tanbeta = get_column_names(tanbeta_names, df)
    m12squared = get_column_names(m12squared_names, df)
    vs = get_column_names(vs_names, df)

    frontcol = [
        Type,
        Lambda1,
        Lambda2,
        Lambda3,
        Lambda4,
        Lambda5,
        Lambda6,
        Lambda7,
        Lambda8,
        vs,
        tanbeta,
        m12squared,
    ]

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