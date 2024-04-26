#!/bin/bash

touch NCI-60_results

awk -F'\t' '$1 ~ /786-0/ && $3 < 0.05 {print $1, $3, $4, $7, $8, $9}' NCI-60_Cancer_Cell_Lines_table.txt >> NCI-60_results
awk -F'\t' '$1 ~ /A498/ && $3 < 0.05 {print $1, $3, $4, $7, $8, $9}' NCI-60_Cancer_Cell_Lines_table.txt >> NCI-60_results
awk -F'\t' '$1 ~ /ACHN/ && $3 < 0.05 {print $1, $3, $4, $7, $8, $9}' NCI-60_Cancer_Cell_Lines_table.txt >> NCI-60_results
awk -F'\t' '$1 ~ /CAKI-1/ && $3 < 0.05 {print $1, $3, $4, $7, $8, $9}' NCI-60_Cancer_Cell_Lines_table.txt >> NCI-60_results
awk -F'\t' '$1 ~ /RXF 393/ && $3 < 0.05 {print $1, $3, $4, $7, $8, $9}' NCI-60_Cancer_Cell_Lines_table.txt >> NCI-60_results
awk -F'\t' '$1 ~ /SN12C/ && $3 < 0.05 {print $1, $3, $4, $7, $8, $9}' NCI-60_Cancer_Cell_Lines_table.txt >> NCI-60_results
awk -F'\t' '$1 ~ /TK-10/ && $3 < 0.05 {print $1, $3, $4, $7, $8, $9}' NCI-60_Cancer_Cell_Lines_table.txt >> NCI-60_results
awk -F'\t' '$1 ~ /UO-31/ && $3 < 0.05 {print $1, $3, $4, $7, $8, $9}' NCI-60_Cancer_Cell_Lines_table.txt >> NCI-60_results
awk -F'\t' '$1 ~ /RXF-631/ && $3 < 0.05 {print $1, $3, $4, $7, $8, $9}' NCI-60_Cancer_Cell_Lines_table.txt >> NCI-60_results
awk -F'\t' '$1 ~ /SN12K1/ && $3 < 0.05 {print $1, $3, $4, $7, $8, $9}' NCI-60_Cancer_Cell_Lines_table.txt >> NCI-60_results