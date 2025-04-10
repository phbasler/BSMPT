#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

import pandas as pd
import sys


def convert(InputFile, OutputFile):
    print(f"Reading {InputFile}.")
    print(f"Output is saved to {OutputFile}.")

    with open(InputFile, "r") as file:
        df = pd.read_csv(file, index_col=False, sep="\t")

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
        print(f"Reb1 changed to {Reb1}.")
        NoImb1 = True
    if Ima1 not in df:
        Rea1 = "a1"
        print(f"Rea1 changed to {Rea1}.")
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

    with open(OutputFile, "w") as file:
        df.to_csv(file, index=True, sep="\t")


if __name__ == "__main__":
    convert(sys.argv[1], sys.argv[2])
