# Criticisms of MAGIC in scRNA-seq data

This repository contains code to explore some of the consequences of using MAGIC in scRNA-seq data.
The `report` subdirectory contains LaTeX files for the report.
The `scripts` subdirectory contains R code to reproduce the analyses described in the report:

- `nostruct_sim.R`, to examine the effect of using MAGIC in a simulation without any structure.
- `structloss_sim.R`, involving a simulation containing primary and secondary separations.
- `fakede_sim.R`, involving DE analysis between two well-separated subpopulations.
- `run_ercc.R`, to examine the effect of using MAGIC on the ERCC data set from 10X Genomics.
- `run_293t.R`, to examine the effect of using MAGIC on the 293T data set from 10X Genomics.
- `run_tb.R`, to test behaviour when combining B and T cells.
