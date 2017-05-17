# differentiation_speciation

Included in this data package are scripts and files to replicate analyses in:

Harvey MG, Seeholzer GF, Smith BT, Rabosky DL, Cuervo AM, Brumfield RT. In press. Positive association between population genetic differentiation and speciation rates in New World birds. *Proceedings of the National Academy of Sciences*.

Included are:

- 100_complete_trees.zip

A file including 100 trees from the pseudo-posterior of Jetz et al. (2012).

- Hackett_lumped_justgeneticdata.tre
- Hackett_split_justgeneticdata.tre
- Hackett_time_justgeneticdata.tre

Trees pruned to each taxonomy examined in the paper, containing just species represented by genetic data in Jetz et al. (2012).

- DR_calculation_lumped
- DR_calculation_split
- DR_calculation_time

R scripts calculating the DR statistic (Jetz et al. 2012) for trees corresponding to the three taxonomies examined in the paper.

- Hackett_lumped_eventsample.rda
- Hackett_split_eventsample.rda
- Hackett_time_eventsample.rda

Rdata files containing event samples from BAMM analyses of the Jetz et al. (2012) trees pruned to reflect each taxonomy.

- lumped_data_final.txt
- split_data_final.txt

Data on population differentiation, species age, and ecological/environmental traits in each species from Smith et al. (2017). 

- taxonomy_key.txt

A key relating species names from our study to those of Jetz et al. (2012).

- Main_text_analyses.R
- Supplementary_analyses.R

Scripts to conduct the analyses used in the main text and supplementary text, respectively.

References:

Jetz W, Thomas GH, Joy JB, Hartmann K, Mooers AO (2012) The global diversity of birds in space and time. Nature 491:444–448.

Smith BT, Seeholzer GF, Harvey MG, Cuervo AM, Brumfield RT (2017) A latitudinal phylogeographic diversity gradient in birds. PLoS Biol 15:e2001073.