# RNA based silent sites

## Overview
The code provided in this repository serves the purpose of identifying and ranking transcriptionally silent sites using RNA sequencing data. The code itself uses Python and Python libraries. Mapping is being performed using BWA mem2 and Samtools.

## Features
- identification of silent sites using RNA sequencing data and reference genomes incl. annotation
- creates convenient output table and report figure
- paired-end data is supported

## Important
This code has been developed to conduct analyses as outlined in the publication above and is currently in development. The code might contain bugs and is also not written for speed or for readability. Future releases might improve on those aspects.

## Disclaimer
This software is intended for research purposes only and is provided "as is" without warranty of any kind. The developer(s) and provider(s) of this software make no representations or warranties, express or implied, regarding the use or performance of this software.

You assume full responsibility and risk for the use of this software. The developer(s) and provider(s) shall not be liable for any direct, indirect, incidental, special, exemplary, or consequential damages arising out of the use of this software.

This software is not intended to be used for critical, medical, or life-saving purposes. It should not be relied upon as the sole basis for decision-making.

By using this software, you acknowledge and agree to this disclaimer. If you do not agree to these terms, you should refrain from using the software.

Last updated: 2024-01-08

## Installation
### Pre-requisites (versions are suggestions and others may work)
- Python 3.9 >= (https://www.python.org/downloads/release/python-390/)
- BWA mem2 >= 2.2 (https://github.com/lh3/bwa)
- samtools >= 1.15 (https://github.com/samtools/samtools)
- The following additional Python libraries are required and can be installed through PIP (see below):
  - Pysam >= 0.22.0
  - Numpy
  - Matplotlib

### bash Setup  
``````
git clone https://github.com/MPUSP/lautenschlaeger_silent_sites.git  
cd lautenschlaeger_silent_sites  
``````
It is recommended to use Python's virtual environments to run the code and to keep your system clean:  
``````
python3 -m venv ./venv  
./venv/bin/pip3 install -r requirements.txt
``````

## Usage

### Data preparation
You will have to prepare input files in tab-separated table format in order to run the analysis.  
Two files are required:
- Read input table:
  - 2-3 columns, tab separated
  - 1st column: path to reference fasta (can be gzip) for this specific NGS dataset (might be repeated with multiple datasets using the same reference)
  - 2nd (and 3rd) column: (paired-end) read file(s) in fastq format (can be gzip)
- Annotation input table:
  - 2-3 columns, tab separated
  - 1st column: path to reference fasta (can be gzip), needs to be repeated when providing multiple information for the same reference
  - either only 2nd column with path to gff file (can be gzip)
  - or column 2 and 3:
    - 2nd column with category of additional annotation, i.e. "prophage"
    - 3rd column with coordinates separated by ":", i.e. 10345:12345, inclusive, 0-based
- each fasta reference should contain only 1 sequence 


### Config file
The tool will try to find executables for bwa and samtools in your path variables. If this is not the case, you can provide a config file with the paths to the executables. With the config file, you will also be able to specify the version of the programs to use.
The config file is a tab-separated text file with 2 columns and currently 2 rows to import program paths into the analysis tool.
  - bwa &nbsp;&nbsp;&nbsp;&nbsp; absolute-path-to-bwa-mem2-executable
  - samtools &nbsp;&nbsp;&nbsp;&nbsp; absolute-path-to-samtools-executable

### Parameters:
| Parameter | Description [default] |
| --- | --- |
| -h | help |
| -c | config file [] (see above) |
| -n | number of cores [1] |
| -i | input table [] (see above) |
| -t | temporary folder [./temp] |
| -o | output folder [./output] |
| -r | minimum region size [100] |
| -p | assumed promoter region [250] |
| --anno | annotation file table [] (see above) |
| --coverage_threshold | normalized coverage threshold [0.02] |
| --verbose | verbose mode (more talkative) |


### Execution:
``````
<path-to-your-virtual-environment>/venv/bin/python3 <path-to-the-code>/silent_sites_analysis.py <parameters to run>
``````

### Example:
We are assuming, that all input files are available at those relative paths.
``````
venv/bin/python3 ./silent_sites_analysis.py  -n 10 -i test/input_table.txt --anno test/anno_table.txt --verbose
``````


## Copyright and License
This software has been created by Knut Finstermeier, Max Planck Unit for the Science of Pathogens, Berlin, Germany, member of the Max Planck Society, https://www.mpg.de/en. This software is released under the GNU General Public License v3.0. See the [LICENSE](LICENSE) file for details. Please cite our software if you use it in your work (see below).

## Release:

Version 1.0, 2024-01-08

## Acknowledgments
This software was written in the context of the referenced publication by Nina Lautenschläger et al. (2024). The authors would like to thank the Max Planck Unit for the Science of Pathogens for support. 

## Citations
Please cite the following publication if you use this software in your work:

TBD

## Contact
Please use software-dev@mpusp.mpg.de for any questions or comments.