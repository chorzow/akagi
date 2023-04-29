# akagi
Name is derived from my favourite anime character but if you really want some meaning in this acronym, it can stand for Archaeal (proKAryotic) Genomic Instruments. Or something else.

## DISCLAIMER
As one can notice from the package name, it is suitable for prokaryotes. Eukaryotic data processing is still WIP and, frankly, I don't know if I really need to implement this modules because there is no urge for that.

## Brief introduction
This Python package simplifies work with Hi-C data by providing numerous useful functions for:
 * launching the most popular pipelines to get contact matrices from raw Hi-C reads 
 * parsing data after calculations are finished (descriptive statistics, contact matrices, etc)
 * exploratory data analysis and visualization

## Structure

For now, akagi consists of three main parts:
 1. `pipelines` contains scripts for launching popular Hi-C pipelines: juicer, distiller and HiC-Pro. This module is run from the command line.
 2. `afterwork` launches tools for downstream analysis such as compartmentalization analysis, TAD and loop calling, RedC signal annotation. This module is run from the command line.
 3. `imagineer` plots the data in suitable format. This module was designed to be used in Jupyter Notebook.
