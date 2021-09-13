# eXNVerify

eXNVerify (Exon and SNV verifier) are the Python tools for extraction and verification of genome sequence fragments coverage quality and presents the results of analysis in an intuitive way for genetic diagnostician. Two executables from this repository takes mosdepth- or bedtools-generated BED file as the whole genome/exome sequence coverage and are able to:
1. (geneCoverage.py) generates detailed verfication of pathogenic germline and somatic single nucletide variants for chosen gene(s)
2. (snvScore.py) analyses the whole genome sequence coverage and evalute all pathogenic germline and SNV coverage quality

Both tools require BED file with the general coverage of WGS/WES sample. In the actual project mosdepth as a fast tool for BAM file analysis was utilized. Detailed description of mosdepth can be found in <link> 


