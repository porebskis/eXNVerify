# eXNVerify

eXNVerify (Exon and SNV verification) are the Python tools for extraction and verification of genome sequence fragments coverage quality and presents the results of analysis in an intuitive way for genetic diagnostician. Two executables from this repository takes mosdepth- or bedtools-generated BED file as the whole genome/exome sequence coverage and are able to:
1. (geneCoverage.py) generates detailed verfication of pathogenic germline and somatic single nucletide variants for chosen gene(s)
2. (snvScore.py) analyses the whole genome sequence coverage and evalute all pathogenic germline and SNV coverage quality

Both tools require BED file with the general coverage of WGS/WES sample. In the actual project mosdepth as a fast tool for BAM file analysis was utilized. Detailed description of mosdepth can be found in TODO.

## Installation

eXNVerify is implemented with Python 3.8 and utilizes popular libraries as `numpy`, `pandas`, and `matplotlib`. It is prepared for Unix operating systems and it is available to pull from DockerHub and run as a standalone Docker container:

`docker pull porebskis/exnverify:0.89b`

## Usage

Repository contains two standalone Python executable codes:

1. geneCoverage.py
```
geneCoverage [-h]
                    SampleBED RefExomeBED SNVGermlineTXT SNVSomaticTXT
                    Threshold GeneName_s [GeneName_s ...]

positional arguments:
  SampleBED       Path to the mosdepth per-base BED output
  RefExomeBED     Path to the all exons BED file
  SNVGermlineTXT  Path to Clivar-generated table with pathogenic germline SNVs
  SNVSomaticTXT   Path to Clivar-generated table with pathogenic somatic SNVs
  Threshold       Coverage quality threshold
  GeneName_s      Gene name(s)

optional arguments:
  -h, --help      show this help message and exit
```

Exemplar docker run coomand with geneCoverage.py execution as follows:
```
docker run -it --rm -v ~/hostpath/:/input -v ~/hostpath/:/output porebskis/exnverify:0.89b ./ExnVerify.py input/SampleBED input/RefExomeBED input/SNVGermlineTXT input/SNVSomaticTXT Threshold GeneName_s
```

The main output of the ``geneCoverage.py`` is the figure with chromosome region where all exones and SNV position occurs for input gene(s). Below exemplar coverage diagram for BRCA1 gene on suitable fragment of chromosome 17:




3. snvScore.py
```
SNVScore [-h] SampleBED SNVGermlineTXT SNVSomaticTXT [Threshold]

positional arguments:
  SampleBED       Path to the mosdepth per-base BED output
  SNVGermlineTXT  Path to Clivar-generated table with pathogenic germline SNVs
  SNVSomaticTXT   Path to Clivar-generated table with pathogenic somatic SNVs
  Threshold       SNV coverage quality threshold (optional, positive)

optional arguments:
  -h, --help      show this help message and exit
```

