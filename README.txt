DEPENDENCIES:

The segmentation script includes plotting support if matplotlib is installed.
The block bootstrap code can plot the test distribution histogram and qqplots
if rpy is installed.

INSTALL:

Jsut enter the directory and run the script, or move the directory to 
site-packages to use as a module.

TEST BLOCK BOOTSTRAP:

To test the block bootstrap code, run

./block_bootstrap.py \
 -1 test_data/conserved_sequences.bed -2 test_data/ENCODE_annotations.bed \
 -d test_data/ENCODE_region_lengths.txt -r 0.041 -n 20 -v

For a detailed explanation of the command line options, run:

./block_bootstrap.py --help

TEST SEGMENTATION:

NOTE: Segmentation can only decrease p-values.  If you run the block bootstrap, 
are interested only in testing, and find a very small p-value, there is probably 
not any reason to go back and segment the data.  That p-value will only get 
smaller after segmentation.  However, if you are building confidence intervals, 
it may still be useful.

To test the segmentation code, first:

-generate test data

./segmentation.py -g

-test the segmentation on the generated data

./segmentation.py -1 sim_output_1.bed -2 sim_output_2.bed \ 
    -d sim_lengths.txt -m 5000 -s 4 -v -p

For a detailed explanation of the command line options, run:

./segmentation.py --help

INPUT FILE FORMATS:

Feature Input Files:

The feature input files must contain new line delimited regions of the form:

chrom ( whitespace ) chromStart ( whitespace ) chromEnd ( whitespace ) ( any other fields )\n

where the 3 necessary fields are:

field 	    type 	description

chrom 	    string 	Name of the chromosome
chromStart 	int 	Starting position on the chromosome
chromEnd 	int 	Ending position on the chromosome 

This is compatible with all BED(>=3)+X formats. More specifically, this is 
compatible with the following ENCODE standard file formats, as described at

http://genomewiki.cse.ucsc.edu/EncodeDCC/index.php/File_Formats

narrowPeak: Narrow (or Point-Source) Peaks Format 
broadPeak: Broad Peaks (or Regions) Format 
gappedPeak: Gapped Peaks (or Regions) Format 
tagAlign: Tag Alignment Format 
pairedTagAlign: Tag Alignment Format for Paired Reads 
NRE Bed6 Format 
BiP Bed8 Format 

Region Length Input Files:

The region lengths input file should consist solely of lines of the form

chrom ( whitespace ) chromStart ( whitespace ) chromEnd \n

where the 3 necessary fields are:

field 	    type 	description

chrom 	    string 	Name of the chromosome
chromStart 	int 	Starting position on the chromosome
chromEnd 	int 	Ending position on the chromosome 

The file format is described at 

http://genomewiki.cse.ucsc.edu/EncodeDCC/index.php/File_Formats

under the heading genomicCoverageFile.

CHANGELOG:

0.5.1 - initial release

marginal and conditional basepair overlap
single file segmentation

0.6.1

added multiple regions segmentation

0.7.1

changed the lengths file format
optimized segmentation for large, binary feature regions
several bug fixes ( see svn logs for more details )

0.7.2

fixed a bug related to 64-bit platform compilation ( Thanks to Michiel de Hoon for the fix )
fixed a bug with specific types of lengths files ( Thanks to Ian Durham for the report )
    sometimes manifest itself as a ZeroDivision Exception,
    sometimes as a negative region length error.
fixed a segmentation file write output bug ( thanks to Mikhail Spivakov for the report )

0.7.3

fixed a mac specific compile bug ( Thanks to Ian Durham for the report )

0.7.4

fixed a python 2.3 64-bit compile bug ( Thanks to Ian Durham for the report )
fixed a python 2.3 compatability bug ( Thanks to Ian Durham for the report )

0.8.0

Rewrite. Added new tests for continuous data, made it easier to add new tests,
and removed some superfluous code.


