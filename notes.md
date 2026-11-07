# My project notes

## Monday 2 Nov: getting the data

- I downloaded public RNA-seq data from GEO, dataset GSE164416 (Wigger et al., Nature Metabolism 2021).
- The samples are islets (the insulin-making part of the pancreas) from 133 patients who had pancreas surgery.
- The data has three parts:
  - the count table: 39,376 genes x 133 samples. Each number = how many reads came from that gene in that sample.
  - the sample information: which group each sample is in (ND, T2D, etc.).
  - the gene annotation: turns gene ID numbers into gene names, e.g. 3630 = INS (insulin).
- Problem I fixed: NCBI sent a "robot check" web page instead of the annotation file. I noticed because the file was 21 KB instead of 4.2 MB, and downloaded it in my browser instead.
- Sanity check: insulin had 440,000 to 1.5 million reads per sample, so the samples really are islets.
- Groups in the full dataset: ND 18, T2D 39, IGT 38, T3cD 35, IFG 3.
- GEO has no sex, age or BMI for these patients. It only has diabetes status, library prep date, and the authors' filter flag.

## Tuesday 3 Nov: preparing the data (script 02)

What each part of script 02 does:

1. Keeps only the two groups I'm comparing: 18 ND and 39 T2D (57 samples). The other 76 are excluded.
2. Matches each sample's information to the right column in the count table, using the GSM sample IDs.
3. Adds the library prep date as "batch", so the model can account for samples prepared on different days.
4. Checks the counts are whole, non-negative numbers, and adds up each sample's total reads (library size).
5. Removes genes with very few reads, because they can't show a reliable difference. 19,396 of 39,376 genes kept.
6. Loads the gene names, and saves everything for the next step.

Numbers from the log:

- Reads per sample (library size): 4.4 to 39.5 million, median 22.2 million.
- The 4.4 million sample is much lower than the rest. I'll check it in the quality step, not delete it.

The authors' filter flag (in_ins_filtered_data_subset):

- TRUE = the authors kept the sample, FALSE = they left it out. 18 of my 57 are FALSE (3 ND, 15 T2D).
- FALSE samples have about half the insulin (median 2.7% of reads vs 5.4%), so the flag probably removes samples with little islet tissue.
- But T2D patients can lose insulin-making cells, so filtering on insulin might remove the most diabetic samples.
- My plan: main analysis on all 57 samples, then repeat it on only the TRUE samples to check the results hold up.

## Wednesday 4 Nov: quality checks (script 03)

What script 03 does: packs the data into a DESeq2 object, plots reads per sample, puts all genes on a fair scale (VST), draws a PCA plot and a sample-similarity heatmap. It only looks; nothing is removed.

What I found:

- Library sizes: 4 samples have under a quarter of the typical reads (2 ND, 2 T2D). Because they're in both groups, low depth won't look like a diabetes difference. DESeq2 adjusts for depth, so I kept them.
- PCA: a map where each dot is a sample and close dots have similar gene activity. PC1 is the summary score that captures the biggest differences (45% here), PC2 the next biggest (14.9%).
- ND and T2D don't separate on the PCA, and T2D samples are more spread out. This matches the original paper, which found T2D islets more varied than healthy ones.
- Prep dates don't separate either, so there's no strong batch effect. I kept batch in the model anyway, to be safe.
- No sample is far enough from the rest (more than 3 SD) to be flagged as an outlier.
- What drives PC1? Not diagnosis and not batch. Insulin content is part of it: PC1 and the share of reads from insulin have a correlation of -0.35 (moderate), so samples further right tend to have less insulin.
- The two samples furthest right (GSM5009304 and GSM5009317, both T2D) have very little insulin (about 0.5% of reads, against a typical 5%) but normal read counts. They're the two odd samples in the heatmap too. Either the laser cut caught little islet tissue, or these patients lost many beta cells. I can't tell which, so I kept them and will check that my top results aren't driven by them.

## Friday 6 Nov: differential expression (script 04)

What script 04 does: for each gene, tests whether it is consistently higher or lower in T2D than in ND, while accounting for batch.

How it works:

1. Corrects for library size with "size factors" (mine ranged from 0.25 to 2.05; the 0.25 sample is the low-read one, scaled up about 4 times).
2. Estimates how noisy each gene naturally is between people, so only differences bigger than that count.
3. Tests each gene: a negative binomial model (a count-friendly regression, count ~ batch + condition) with a Wald test for the T2D vs ND effect.
4. Corrects for testing 19,396 genes with the Benjamini-Hochberg method (the adjusted p-value, padj).

Key terms:

- LFC (log2 fold change): how much higher or lower a gene is in T2D. +1 = twice as high, -1 = half, 0 = no change. Positive = higher in T2D, negative = lower.
- FDR (false discovery rate): among the genes I call significant (padj < 0.05), about 5% are expected to be false alarms.

Results:

- 19,396 genes tested.
- 2,954 significant at FDR < 0.05: 1,638 higher in T2D, 1,316 lower in T2D.
- 401 of these change at least 2-fold.
- 770 couldn't be tested properly: 753 had too few reads, 17 had one extreme sample.
- Even though ND and T2D didn't separate on the PCA, thousands of genes shift slightly and consistently. That's why the test goes gene by gene.

Top 5 genes:

| Rank | Gene | Higher or lower in T2D | What it does |
|---|---|---|---|
| 1 | DKK3 | Higher | Blocks Wnt signalling. Known as a tumour suppressor in cancer; in islets the Wnt pathway helps control beta-cell growth and function. |
| 2 | ARG2 | Lower | Enzyme that breaks down L-arginine into ornithine and urea; linked to nitric oxide production. |
| 3 | SLC2A2 | Lower | Glucose transporter 2 (GLUT2), which beta cells use to take in glucose and sense blood sugar. Less of it could mean beta cells sense glucose worse, or that there are fewer beta cells. |
| 4 | LOC105376812 | Higher | Uncharacterised gene, no known function. LFC 5.6 (about 48 times higher) is suspiciously large; check whether only a few samples drive it. |
| 5 | RNF6 | Lower | E3 ubiquitin ligase: tags proteins for breakdown; involved in cell development and cancer progression. |

Pattern and comparison with published work:

- ARG2 (my #2) is one of the 19 genes in a published T2D islet signature (Solimena et al., Diabetologia 2018), where it helped insulin secretion. That study also measured SLC2A2.
- That study came from the same research group and similar patients, so this shows consistency, not fully independent proof.
- Two of my top three (ARG2, SLC2A2) are beta-cell function genes and both are lower in T2D, which fits failing beta cells. It could also partly reflect fewer beta cells in the tissue.
- LOC105376812 is about 48 times higher in T2D, which is suspiciously large. Check on Saturday whether only a few samples drive it.

## Saturday 7 Nov: P-value histogram, MA plot, Volcano plot, and Heatmap

Two patterns separated T2D from non-diabetic islets. Genes involved in beta-cell function, such as ARG2 and SLC2A2, were modestly but consistent lower in T2D across patients.
A second group, including antibody(immunoglobulin) and inflammation-related genes, was much higher in a subset of T2d patients.
This suggests that immune cells or inflammation in some samples rather than a change shared by every patient.

Comparison with published work: I compared my results with a published 19-gene T2D islet signature (Solimena et al., Diabetologia 2018).
14 of the 19 genes were significant in my data, all in the same direction, including ARG2, SLC2A2, TMEM37 and PPP1R1A. Of the other five, four went the same direction without reaching significance (FFAR4 just missed, padj 0.059), and one (FAM102B) was filtered out for low counts.
None went the opposite way. The two studies share a research group and probably some patients, so this shows consistency rather than independent confirmation.
