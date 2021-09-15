# eXNVerify

eXNVerify (Exon and SNV verification) includes Python-based tools for extraction of genome sequence fragments and verification of coverage quality. Tools wrapped into ``Docker container`` present the results of analysis in an intuitive way for genetic diagnostician. Two executables take BED file as the whole genome/exome sequence coverage and are able to:
1. (geneCoverage.py) perform detailed verfication of pathogenic germline and somatic single nucletide variants (SNV) for chosen gene(s),
2. (snvScore.py) analyse the whole genome sequence coverage and evaluate all pathogenic germline and somatic SNV coverage quality.

Both tools require decompressed BED file with the general coverage of WGS/WES sample. In the actual implementation output of the ``mosdepth`` as a fast tool for BAM file analysis. Output per-base BED file is utilized. Detailed description of mosdepth can be found in TODO.

## Installation

eXNVerify is implemented with Python 3.8 and utilizes popular libraries as `numpy`, `pandas`, and `matplotlib`. It is prepared for Unix operating systems and it is available to pull from DockerHub and run as a standalone Docker container:

`docker pull porebskis/exnverify:0.89b`

## Usage

Repository contains two standalone Python executable codes:

### 1. geneCoverage.py
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

Exemplar ``docker run`` command with geneCoverage.py execution as follows:
```
docker run -it --rm -v ~/hostpath/:/input -v ~/hostpath/:/output porebskis/exnverify:0.89b ./geneCoverage.py \
           input/SampleBED input/RefExomeBED input/SNVGermlineTXT input/SNVSomaticTXT Threshold GeneName_s
```
User may decide what kind of source information would be processed by ``geneCoverage`` ie. exome reference BED may contain position of exones related to all genes but it is also possible to choose desired subset depending on the need. User may also choose which SNV should be the inputs of the analysis. For the purpose of this project development, two Clinvar-generated tables are prepared separately for germline and somatic pathogenic SNVs. These tables can be easily generated from the ClinVar. In the repository, user may find examples of these tables: they contain all pathogenic SNV for all genes, however, they are limited (using Clinvar filters) to these records which are approved by expert panels or submitted by multiple sources to Clinvar. Thus, two SNV tables saved as txt files contain about 15k SNVs and user may feel free to utilize them in their samples analysis. Optional argument of ``geneCoverage`` is the user-defined threshold value which indicates user-desired coverage quality of the analysed sample. 

The main output of the ``geneCoverage.py`` is the PDF figure with chromosome region where all exones and SNV position occurs for given gene(s). The code execution prepares at least one PDF figure for each input gene. If exones related to the input gene are located in more than one chromosomee (e.g. PSRR2), suitable number of figures is generated. If the input gene is not listed in ``RefExomeBED`` file, analysis is not performed. 

In this way, ``geneCoverage.py`` allows verifying the gene coverage in detail and support diagnostician to evaluate the wgs/wes sample for the purposes of further analysis and diagnosis process.

### 2. snvScore.py
```
snvScore [-h] SampleBED SNVGermlineTXT SNVSomaticTXT [Threshold]

positional arguments:
  SampleBED       Path to the mosdepth per-base BED output
  SNVGermlineTXT  Path to Clivar-generated table with pathogenic germline SNVs
  SNVSomaticTXT   Path to Clivar-generated table with pathogenic somatic SNVs
  Threshold       SNV coverage quality threshold (optional, positive)

optional arguments:
  -h, --help      show this help message and exit
```

The ``docker run`` command with ``snvScore.py``execution is as follows:
```
docker run -it --rm -v ~/hostpath/:/input -v ~/hostpath/:/output porebskis/exnverify:0.89b ./snvScore.py input/SampleBED input/SNVGermlineTXT input/SNVSomaticTXT Threshold
```
Similar to the ``geneCoverage``, ``snvScore`` requires Clinvar-generated tables with recors of pathogenic germline and somatic SNVs. User may generate their source tables using Clinvar search tool and filter or utilizes included TXT files in this eXNVerify repository. Optional argument of ``snvScore`` is the user-defined threshold value which indicates user-desired coverage quality of the analysed sample. 

The output of ``snvScore.py`` is the report TXT file with coverage information about all pathogenic single nucleotide variants (germline and somatic) in the input sample. All SNVs can be taken from Clinvar repository, user may choose the SNVs from Clinvar, the authors prepared the input table as the collection of all variants that are related to pathogenic SNV. 

## Example outputs

### 1. geneCoverage

Below exemplar coverage diagram for BRCA1 gene in HG003 sample from pacbio-hifi (downloaded from Google Storage):

![BRCA1 coverage](/fig/BRCA1.chr17.HG003.pacbio-hifi.21x.haplotag.grch38.bam.per-base.bed.png)

Additionally, as summary log:
```
94% of all pathogenic germline SNVs and 98% of all pathogenic somatic SNVs are covered above threshold (15)
```

### 2. snvScore
Exemplar coverage report as the results of ``snvScore.py`` HG003 sample from pacbio-hifi (downloaded from Google Storage) is as follows:
```
SNV coverage report - HG003.pacbio-hifi.21x.haplotag.grch38.bam.per-base.bed

81% of all pathogenic germline SNVs and 91% of all pathogenic somatic SNVs are covered above threshold (15)

Whole genome coverage:
median mean  std 1st quartile 3rd quartile  min     max
    22   32  190           18           26    0   24444

Pathogenic (G - germline, S - somatic) SNV coverage:
('count' is the number of variants in a given region)

region count(G)  median(G) std(G) min(G) max(G) count(S)  median(S) std(S) min(S) max(S)
----------------------------------------------------------------------------------------
   ALL    14157         20      6      1     43      310         20      5      9     36
----------------------------------------------------------------------------------------
  chr1     1182         21      5      8     37       10         24      3     17     29
  chr2     1118         22      6      5     39       33         20      5     12     33
  chr3      774         20      5      8     37       40         21      6     10     33
  chr4      270         20      7      7     41        6         16      2     15     19
  chr5      538         21      5      5     36       10         18      4     15     29
  chr6      511         22      6      4     43        0          0      0      0      0
  chr7      717         22      5      6     43       13         26      4     21     34
  chr8      335         23      6      8     37        0          0      0      0      0
  chr9      479         20      6      2     37        1         27      0     27     27
 chr10      362         20      6      9     36       32         20      7      9     36
 chr11      966         19      6      6     36        9         14      6     11     27
 chr12      709         22      5     10     38       12         23      3     16     25
 chr13      910         17      5      9     37       26         22      5     11     29
 chr14      324         23      4      9     34        6         21      4     13     27
 chr15      533         21      5      8     35        4         17      1     17     20
 chr16      643         19      6      3     40        1         17      0     17     17
 chr17     1561         20      4      8     37       88         21      2     14     24
 chr18      154         20      6     11     35        4         25      2     25     30
 chr19      720         18      4      7     27       10         12      4      9     20
 chr20      146         18      5      8     35        0          0      0      0      0
 chr21      174         14      6      2     28        1         22      0     22     22
 chr22      181         16      5      1     29        2         18      6     12     24
  chrX      850         10      4      3     24        2         14      4      9     18
  chrY        0          0      0      0      0        0          0      0      0      0
----------------------------------------------------------------------------------------
```
The 3rd line of report file is also a summary log send to the standard output when ``snvScore`` is finished with success.
