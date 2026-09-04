# Input formats

## MGF

Each spectrum must contain `BEGIN IONS`, `END IONS`, `SCANS`, `PEPMASS`, and
either `RTINMINUTES` or `RTINSECONDS`. `SCANS` must match the feature identifier
used in the abundance table.

## Abundance table

CSV or TSV with the following exact columns:

```text
Alignment ID,Average Rt(min),Average Mz,A,B,C,D,E
```

Rows are aligned precursor features. Columns A-E are non-negative feature
abundances. Replicates should be summarized consistently before this step.

## Activity-response table

Long CSV or TSV:

```text
Group,Response
Model,13.79
Control,2.50
A,14.54
...
```

Groups `Model`, `Control`, and `A`-`E` are required. Multiple replicate rows per
group are supported. Alternative column names can be passed with
`--activity-group-column` and `--activity-response-column`.

## Existing edge table

Required columns:

```text
Source_Scan,Target_Scan,Matched_Peaks,Modified_Cosine_Score
```

Precursor columns are optional because they can be recovered from the MGF.

## Candidate structures

Long CSV with one row per candidate and these exact columns:

```text
Alignment_ID,Name,Formula,SMILES,Ontology,Score,Rank
```

Only parsable structures whose formulas contain C/H/O and optionally N are used.
