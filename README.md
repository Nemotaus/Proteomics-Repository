# Proteomics-Repository
The R Code for the Proteomics Data in Figure 2 of Repairing Atp10D in C57Bl/6J mice restores protein expression but does not mitigate metabolic stress from high fat diet.

This code was given to me by J. Hermanson for my manuscript Repairing Atp10D in C57Bl/6J mice restores protein expression but does not mitigate metabolic stress from high fat diet. However, I made significant edits.

This code was used to merge two excel data sheets from two independent mass spectrometry runs. Supplemental File 1 is the merged Excel sheet. Though be warned my dearest scallywags, it is possible the "spaces", +/+, and -/-'s could affect R's ability to identify columns. I editted this Excel sheet so it made sense for anyone to be able to see the individual mice and their genotypes easily, so that they would understand that I used an N of 4. My hope was that this would allow anyone to scour the data in excel to their own enjoyment, although the data is boring in my honest opinion. i would highly recommend editting excel sheet to replace all spaces in column headers with "_" and all "/" with "". the final names should read something like "3291_ATP10D--".

Supplemental File 1 is daf1 in the code therefore it is post-HeLa cell removal. In my full western blots in supplemental File S3 I have HeLa cells that overexpressed murine ATP10D to help find high abundance peptides. I will attach the unmerged excel sheets to this repository.

This R script contains all math and code used to analyze this data set as best as possible. I did export some graphs and data as SVG to adjust the final aesthetics of these figures to my liking in ways that were easier than coding them in R. Also, I was too lazy to learn how to code graphs with that high specificity.

With Supplemental File 1 chosen as your daf1 by the file.choose() command, you should be able to start the code at line 52.

Do not forget to attach dependencies as needed.
