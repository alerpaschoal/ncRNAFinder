![Badge em Desenvolvimento](https://img.shields.io/badge/build-passing-green)
![Python Version](https://img.shields.io/badge/python-3.10-blue.svg)
# ncRNAFinder
ncRNAFinder is an automatic and scalable system for large-scale data annotation analysis of ncRNAs which use both sequence and structural search strategy for ncRNA annotation.

## Install
To use the ncRNAFinder, it is necessary to install some dependencies and databases. First, the necessary tools are BLAST (version 2.15.0) and INFERNAL (1.1.5). Second, the databases needed are RNAcentral (version 24) and Rfam (version 14.10). Lastly, the Python libraries required are biopython, matplotlib, matplotlib_venn, numpy, pandas, and pip. To install them, simply navigate to the folder containing all the codes and use our automatic installation:

IMPORTANT: The download of the RNAcentral database takes some time, almost 7 hours. 
~~~
sudo sh install_ncRNAFinder.sh
~~~

After using our automated installation, our tool will already be in the PATH, allowing you to use it in any folder.

## Usage
To execute the tool, simply use the following command:
~~~
ncRNAFinder.py -f <input_file> -o <output_name>
~~~

### Mandatory parameters:
~~~
-f|input_file <file_name>                                       Input file in FASTA format

-o|output_name <output_name>                                    Output name to save the results
~~~

### Optional parameters:
~~~
-i|pident <integer>                                       Minimun percentage of identity of BLASTn. (Default: 95)

-c|coverage <integer>                                     Minimun percentage of coverage of BLASTn. (Default: 95)

-t|threads <integer>                                      Number of threads. (Default: 1)

-r|BestHit <1|0>                                          Option to filter only the best result between two strands, based on E-value), 1-yes or 0-no. (Default: 1)
~~~

## Reference

## Contact
To report bugs, to ask for help and to give any feedback, please contact Alexandre R. Paschaol (paschoal@utfpr.edu.br) or Vitor Gregorio (vitor-gregorio@hotmail.com).
