## Helper scripts for input data conversion
The scripts in this directory convert data to the tabular-separated
input format required for the BSMPT model implementations.

The column names are adjusted to match the output column names
used by ScannerS v2.0.0 (https://gitlab.com/jonaswittbrodt/ScannerS).
If input is used with custom column names, their names need to be added to the parameter name lists.

The scripts for a model `MODEL` are run with

```bash
python3 prepareData_MODEL.py -in [INPUTFILE] -out [OUTPUTFILE] [OPTIONAL ARGUMENTS]
```

Required arguments are input and output file name, optional
arguments are the column separator of the input file or the
number of the index column. Run the above command with `-h` or
`--help` for a list of all arguments.
