# Proteomics-Repository
The R Code for the Proteomics Data in Figure 2 of Repairing Atp10D in C57Bl/6J mice restores protein expression but does not mitigate metabolic stress from high fat diet.

This code was given to me by J. Hermanson for my manuscript Repairing Atp10D in C57Bl/6J mice restores protein expression but does not mitigate metabolic stress from high fat diet. However, I made significant edits.

This code was used to merge two excel data sheets from two independent mass spectrometry runs. Supplemental File 1 is the merged Excel sheet. However, I editted Supplemental File 1 to be more reader friendly. I attached the original daf1_exported csv to this repository if you want to run the code yourself. I have also uploaded the original DIANN log2 values that are used in both Excel sheet of Supplemental File 3 if you wish to start the code at the merge function.

In the code I write arguments to drop the HeLa column. These arguments are irrelevant to the manuscript and refer to samples not used in the preparation of this manuscript in Supplemental Figure 3. These lines are only needed if you start with the uploaded separate csv files.

This R script contains all math and code used to analyze this data set as best as possible. I did export some graphs and data as SVG to adjust the final aesthetics of these figures to my liking in ways that were easier than coding them in R. Also, I was too lazy to learn how to code graphs with that high specificity.
