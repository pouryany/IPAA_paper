[![DOI](https://zenodo.org/badge/705339111.svg)](https://doi.org/10.5281/zenodo.13895686)

## Overview

This repository provides the code companion for the manuscript entitle:
“Identification of p38 MAPK as a Major Therapeutic Target for
Alzheimer’s Disease based on Integrative Pathway Activity Analysis and
Validation in 3D Human Cellular Models.”

please contact `pnaderiy [at] bidmc [dot] harvard [dot] edu` for any
questions you may have about the codes.

## Citations

## Code instructions

First, you need to access processed RNA-Seq datasets from the ROSMAP,
MSBB, and MAYO cohort via “Synape.org” (**data accessions to be
provided**). These datasets require an approved Data Use Agreement to
protect human subject privacy. You also need to download RNA-Seq
datasets corresponding to the cellular models (**data accession to be
provided**)

The RNA-Seq datasets and associated covariates should be placed into
respective directories under the folder `preprocessed`.

All codes relating to RNA-Seq data processing are placed in the `codes`
directory.

<table>
<colgroup>
<col style="width: 11%" />
<col style="width: 44%" />
<col style="width: 44%" />
</colgroup>
<thead>
<tr class="header">
<th>Folder name</th>
<th>Description</th>
<th>Details</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><code>codes</code></td>
<td>Associated R-script</td>
<td>contains 4 subdirectories that organize codes for cleaning, pathway
dysregulation analysis, plotting, and processing additional
datasets.</td>
</tr>
<tr class="even">
<td><code>figures</code></td>
<td>figure outputs.</td>
<td>Figures are saved in this directory according to their order in the
manuscript</td>
</tr>
<tr class="odd">
<td><code>Output</code></td>
<td>Processed datasets</td>
<td>Dysregulated pathway profiles, pathway activity residuals, and
similarity analysis results</td>
</tr>
<tr class="even">
<td><code>preprocessed</code></td>
<td>Necessary data files</td>
<td>preprocessed RNA-Seq data from cellular models and human datasets.
Background pathways and gene-level information are also provided to
facilitate appropriate codes.</td>
</tr>
<tr class="odd">
<td><code>tables</code></td>
<td>output tables</td>
<td>dysregulated pathways, shared pathways across human cohorts and cell
models, differentially expressed genes.</td>
</tr>
<tr class="even">
<td><code>reduced_pathway_output</code></td>
<td>supplementary experiments</td>
<td>dysregulated pathways produced using a reduced background dataset
for robustness analysis.</td>
</tr>
</tbody>
</table>

## code instructions

<table>
<colgroup>
<col style="width: 11%" />
<col style="width: 44%" />
<col style="width: 44%" />
</colgroup>
<thead>
<tr class="header">
<th>Folder name</th>
<th>file name</th>
<th>description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><code>00_new_msigdb</code></td>
<td><code>01-MSigDBV6CONSERVATIVE.R</code></td>
<td>Supplementary analysis: generating a reduced background pathway
dataset for robustness analysis.</td>
</tr>
<tr class="even">
<td><code>00_new_msigdb</code></td>
<td><code>02_MSigDB_XML.R</code></td>
<td>Generating display friendly names for MSigDB Pathways</td>
</tr>
<tr class="odd">
<td><code>01_cleaning_and_prep</code></td>
<td><code>01-DEG_normalizer.R</code></td>
<td>processing DEG files from multiple human cohorts to provide a
concordant representation</td>
</tr>
<tr class="even">
<td><code>01_cleaning_and_prep</code></td>
<td><code>02-clean_DEG_Tables.R</code></td>
<td>processing DEG files from multiple cell lines to provide a
concordant representation</td>
</tr>
<tr class="odd">
<td><code>02_Dysregulated pathways</code></td>
<td><code>01-FullPipeline_Mayo_ADvsControl.R</code></td>
<td>Dysregulated pathway analysis in the Mayo cohort via PanomiR
package</td>
</tr>
<tr class="even">
<td><code>02_Dysregulated pathways</code></td>
<td><code>02-FullPipeline_MSBB_ADvsControl.R</code></td>
<td>Dysregulated pathway analysis in the MSBB cohort via PanomiR
package</td>
</tr>
<tr class="odd">
<td><code>02_Dysregulated pathways</code></td>
<td><code>03-FullPipeline_ROSMAP_ADvsControl.R</code></td>
<td>Dysregulated pathway analysis in the Mayo cohort via PanomiR
package</td>
</tr>
<tr class="even">
<td><code>02_Dysregulated pathways</code></td>
<td><code>04-FullPipeline_A5vsG2B2.R</code></td>
<td>Dysregulated pathway analysis in the A5 cells via PanomiR
package</td>
</tr>
<tr class="odd">
<td><code>02_Dysregulated pathways</code></td>
<td><code>05-FullPipeline_D4vsG2B2.R</code></td>
<td>Dysregulated pathway analysis in the D4 cells line via PanomiR
package</td>
</tr>
<tr class="even">
<td><code>02_Dysregulated pathways</code></td>
<td><code>06-FullPipeline_H105vsG2B2.R</code></td>
<td>Dysregulated pathway analysis in the h10 cells via PanomiR
package</td>
</tr>
<tr class="odd">
<td><code>02_Dysregulated pathways</code></td>
<td><code>07-FullPipeline_I47FvsG2B2.R</code></td>
<td>Dysregulated pathway analysis in the I47F cells via PanomiR
package</td>
</tr>
<tr class="even">
<td><code>02_Dysregulated pathways</code></td>
<td><code>08-FullPipeline_I45FvsI47F.R</code></td>
<td>Dysregulated pathway analysis in the I45F cells compared to I47F
cells via PanomiR package</td>
</tr>
<tr class="odd">
<td><code>02_Dysregulated pathways</code></td>
<td><code>09-clean_Pathways.R</code></td>
<td>cleaning up dysregulated pathway tables</td>
</tr>
<tr class="even">
<td><code>02_Dysregulated pathways</code></td>
<td><code>10-Gene_based_correlation.R</code></td>
<td>Correlation analysis of gene dysregulation across cell models and
human cohorts</td>
</tr>
<tr class="odd">
<td><code>02_Dysregulated pathways</code></td>
<td><code>11-Pathway_based_correlation2_new.R</code></td>
<td>Correlation analysis of gene dysregulation across cell models and
human cohorts. Heatmap of similarity in Figure 3D</td>
</tr>
<tr class="even">
<td><code>03_plots</code></td>
<td><code>01-assay_similarities_brains.R</code></td>
<td>Venn diagrams, heatmaps, and correlation plots in Figure 2. Tables
representing shared dysregulated pathways. Venn diagram in Figure
4.</td>
</tr>
<tr class="odd">
<td><code>03_plots</code></td>
<td><code>02-region_plotting.R</code></td>
<td>Pathway activity heatmap in Figure 2</td>
</tr>
<tr class="even">
<td><code>03_plots</code></td>
<td><code>03-Maximal_Heatmap.R</code></td>
<td>Pathway activity heatmap in Figure 5</td>
</tr>
<tr class="odd">
<td><code>03_plots</code></td>
<td><code>04-Bar_plots.R</code></td>
<td>Heatmaps of shared dysregulated pathways in Figure 4. Co-expression
network plot</td>
</tr>
<tr class="even">
<td><code>03_plots</code></td>
<td><code>05_PCA_covariates.R</code></td>
<td>PCA-based determination of significant confounding variables in
Figure S01.</td>
</tr>
<tr class="odd">
<td><code>03_plots</code></td>
<td><code>06-gene_correlation_similarities_2.R</code></td>
<td>Correlation analysis of differentially expressed genes, Figure
2F</td>
</tr>
<tr class="even">
<td><code>03_plots</code></td>
<td><code>07-get_ma_plots.R</code></td>
<td>MA plots corresponding to Figure S02</td>
</tr>
<tr class="odd">
<td><code>03_plots</code></td>
<td><code>09-assay_similarities_new.R</code></td>
<td>Bar plots associated with Chi-squared tests in Figure 3C</td>
</tr>
<tr class="even">
<td><code>03_plots</code></td>
<td><code>10_Batch_correction.R</code></td>
<td>PCA plot in Figure S02</td>
</tr>
<tr class="odd">
<td><code>03_plots</code></td>
<td><code>11_get_P38_Genes.R</code></td>
<td>Activity of genes in the P38 MAPK Pathway in Figure 5A</td>
</tr>
<tr class="even">
<td><code>03_plots</code></td>
<td><code>12_Visualize_gene_matrix.R</code></td>
<td>Gene-based heatmap of similarity in Figure 3E</td>
</tr>
<tr class="odd">
<td><code>03_plots</code></td>
<td><code>S01-survival_brains.R</code></td>
<td>P-value distribution analysis in Figure S02</td>
</tr>
<tr class="even">
<td><code>external datasets</code></td>
<td><code>01_external_dataset_clean.R</code></td>
<td>cleaning iPSC gene expression data</td>
</tr>
<tr class="odd">
<td><code>external datasets</code></td>
<td><code>02_external_similarity_sq.R</code></td>
<td>pathway-based similarity analysis Figure S03</td>
</tr>
<tr class="even">
<td><code>external datasets</code></td>
<td><code>03_External_heatmap.R</code></td>
<td>pathway-based similarity analysis Figure S03</td>
</tr>
</tbody>
</table>

Supplementary files that are necessary for reproducing the study are
provided in the `preprocssed/` directory

<table>
<colgroup>
<col style="width: 11%" />
<col style="width: 44%" />
<col style="width: 44%" />
</colgroup>
<thead>
<tr class="header">
<th>File name</th>
<th>Description</th>
<th>Reference/Resource</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><code>geneParameters.tsv</code></td>
<td>GC-content and gene-length</td>
<td>ENSEMBL</td>
</tr>
<tr class="even">
<td><code>MSigDBPathGeneTab.RDS</code></td>
<td>Background Pathways for the PanomiR package</td>
<td>MSigDB and PanomiR Packages</td>
</tr>
<tr class="odd">
<td><code>MSigDBPathGeneTabLite.RDS</code></td>
<td>Reduced overlap Background Pathways MSigDB, for robustness
testing</td>
<td>MSigDB and PanomiR Packages</td>
</tr>
<tr class="even">
<td><code>MsigDB_jaccard.zip</code></td>
<td>Jaccard overlap between MSigDB package, unzip before using</td>
<td>MSigDB</td>
</tr>
<tr class="odd">
<td><code>MsigDB_display_names.csv</code></td>
<td>clean, displayable names from MSigDB</td>
<td>MSigDB</td>
</tr>
<tr class="even">
<td><code>PCxN_MSigDB_withJaccard.RDS</code></td>
<td>Pathway coexpression network of MSigDB database</td>
<td>MSigDB and PCxN</td>
</tr>
</tbody>
</table>
