# De novo mutations in regulatory elements cause developmental disorders
Code for paper published in March, 2017 analysing the role of de novo mutations in severe developmental disorders (part of the deciphering developmental disorders project).

This analysis is broadly separated into three major parts:
1. Population genetics analysis of conserved non-coding elements (uses the dddMAPS library https://github.com/pjshort/dddMAPS).
2. Testing for burden of de novo mutations in non-coding elements and identifying recurrently mutated elements.
3. Estimating general properties of pathogenic variation in conserved non-coding elements (proportion of sites that are pathogenic with high penetrance and genome-wide estimate of contribution of de novo mutations in regulatory elements to severe DD).

## Population Genetics
The Figure1_Element_PopGen.Rmd notebook in analysis_notebooks has the code used to generate all four panels in Figure 1.

## De Novo Burden
The Figure2_DNMburden.Rmd notebook in analysis_notebooks uses the tri-nucleotide mutation rate model described in Samocha et. al, 2014 to test for an enrichment of DNMs across different regulatory element classes.
The Figure3_recurrent_DNMs.Rmd notebook was used to identify recurrently mutated elements, compared the number observed to the expected under the trinucleotide null model, and to test each element individually against a genome-wide significance cutoff.

## Estimating General Properties of Non-Coding Elements
The Figure4_maximum_likelihood_and_genome_estimate.Rmd notebook contains code to assess the likelihood of what we observe (zero elements at genome-wide significance) across a number of parameters that define the pathogenicity 'landscape' of a set of elements as well as code to do the genome-wide extrapolation of DNM burden based on evolutionary conservation score.

This paper has been posted on biorXiv and the code/analyses will continue to be updated and revised.
