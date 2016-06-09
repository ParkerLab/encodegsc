#!/usr/bin/python

import getopt, sys
from base_types import *

def regions_complement(regions):
    comp_regions = base_types.regions()

    for region_name in regions.keys():
        # store an alias for the region with name region_name - just readability
        data = regions[region_name]
        comp_regions[region_name] = base_types.feature_region( data.length, (), region_name )

        # for each feature interval, store the non-region *before* it
        prev_end = -1
        for fi in data.iter_feature_regions():
            if fi.start - prev_end > 1:
                comp_regions[region_name].add(interval(prev_end+1, fi.start-1))
            prev_end = fi.end
         
        # store the non-feature interval after the last feature interval
        if prev_end < data.length - 1:
            comp_regions[region_name].add(interval(prev_end+1, data.length-1))

    return comp_regions

def usage():
  print >> output, """
  Take the complement of a feature set.

  ie, mark every feature a non-feature and every non-feature a feature. 

  REQUIRED ARGUMENTS:
  -i --input : the input file. 
  
  -d --domain : the domain of the data files, the portion of the genome over 
  which the features ( from --input ) are defined.

  The file formats are described in the README under INPUT FILE FORMATS.

  OPTIONAL ARGUMENTS:

  -p --prefix: The prefix for the new output files. Defaults to 'complement_'.

  -v --verbose: Print extra output

  -h --help: Print this help info and exit.
"""

def main():
    try:
        long_args = ["input=", "domain=", "output_file_prefix=", "verbose", "help"]
        opts, args = getopt.getopt(sys.argv[1:], "i:d:o:vh", long_args)
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

    # set the default options
    output_fname_prefix = "complement_"

    for o, a in opts:
        if o in ("-v", "--verbose"):
            global verbose
            verbose = True
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-i", "--input"):
            bed_fp = open(a)
        elif o in ("-d", "--domain"):
            domain_fp = open(a)

    try: 
        assert vars().has_key('bed_fp')
        assert vars().has_key('domain_fp')
    except Exception, inst:
        usage()
        print inst
        sys.exit()
    
    regions = parse_bed_file(bed_fp, domain_fp)
    bed_fp.close()
    domain_fp.close()

    regions_comp = regions_complement(regions)
    # free this memory in case we need it for the file write    
    del regions

    of = open(prefix_filename(bed_fp.name, output_fname_prefix), 'w')
    olf = open(prefix_filename(domain_fp.name, output_fname_prefix), 'w')
    regions_comp.writeBedFile(of, olf)
    of.close()
    olf.close()

if __name__ == "__main__":
    main()

