# Macrophages

This repo covers the steps taken to look for candidate genes in Endoplasmic Reticulum (ER) membrane that can be induced/inhibit in order to promote M1-macrophage polarization.

To see the workflow and follow along, click this page: https://gsanramon.github.io/macrophages/GSE57614.mRNA.data.html

The following files are included in this repo:
1. GSE57614.mRNA.data.rmd
2. GSE57614.mRNA.data.html
3. GSE57614.mRNA.data.R - standalone R script that can be used to integrate with a pipeline 
4. Input files:
  + ER.txt - annotations of ER proteins
  + GSE57614.M0vsM1.24hr.txt - M0 vs M1 (24hrs) RNA data
  + GSE57614.M2vsM1.24hr.txt - M2 vs M1 (24hrs) RNA data
  + Ab_Chris_refined.csv - antibody data
  + phospho_Chris_refined.csv - antibody data
5. Output files in a separate folder
  + common2 - list of genes that are differentially expressed and common to both M0vsM1 and M2vsM1
  + sort.M0vsM1.DEG.ERmem - data of candidate genes from M0 vs M1 file 
  + sort.M2vsM1.DEG.ERmem - data of candidate genes from M2 vs M1 file
