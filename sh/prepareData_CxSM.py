#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
# Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

import pandas as pd
import argparse


def convert(IndexCol, InputFILE, OutputFILE, Seperator):
    df = pd.DataFrame()
    if IndexCol == "False":
        df = pd.read_table(InputFILE, index_col=False, sep=Seperator)
    else:
        df = pd.read_table(InputFILE, index_col=int(IndexCol), sep=Seperator)

    """
    The parameters should have the label of the corresponding parameter.
    """
    v = "v"
    vs = "vs"
    va = "va"
    msq = "msq"
    lamb = "lambda"
    delta2 = "delta2"
    b2 = "b2"
    d2 = "d2"
    Reb1 = "Reb1"
    Imb1 = "Imb1"
    Rea1 = "Rea1"
    Ima1 = "Ima1"

    NoImb1 = False
    NoIma1 = False

    if Imb1 not in df:
        Reb1 = "b1"
        NoImb1 = True
    if Ima1 not in df:
        Rea1 = "a1"
        NoIma1 = True

    frontcol = [v, vs, va, msq, lamb, delta2, b2, d2, Reb1, Imb1, Rea1, Ima1]

    for c in frontcol:
        if c not in df:
            df[c] = 0

    Col = [c for c in frontcol if c in df] + [c for c in df if c not in frontcol]

    df = df[Col]

    if NoImb1:
        df.rename(columns={Reb1: "Reb1"}, inplace=True)
    if NoIma1:
        df.rename(columns={Rea1: "Rea1"}, inplace=True)

    with open(OutputFILE, "w") as file:
        df.to_csv(file, index=True, sep="\t")


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