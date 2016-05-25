# /usr/bin/env python
'''
This script will convert an input VCF file into an output file which can be passed to APOLLOH for the allelic count
information. The script will read through the VCF file and identify heterozygous SNPs. Data for these positions will
then be parsed and written to the output file in the appropriate format.

Dependencies:
    PyVCF >= 0.6.3

@author Andrew Roth
@license GPLv3
'''
import csv
import vcf

def main(args):
    normal_column = args.normal_column
    
    if normal_column == 0:
        tumour_column = 1
    else:
        tumour_column = 0
    
    
    reader = vcf.Reader(open(args.vcf_file))
    
    writer = csv.writer(open(args.out_file, 'w'), delimiter='\t')
    
    for record in reader:
        normal_call = record.samples[normal_column]
        
        if not normal_call.is_het:
            continue
        
        if normal_call.data.GQ < args.min_genotype_quality:
            continue
        
        tumour_call = record.samples[tumour_column]

        # Parse out the required fields
        chrom = record.CHROM
        
        coord = record.POS
        
        ref_base = record.REF
        
        var_base = record.ALT[0]
        
        # Sites with no coverage will be skipped
        try:
            ref_counts = tumour_call.data.AD[0]
        
            var_counts = tumour_call.data.AD[1]
        
        except:
            continue
        
        # Create output row and write it to file
        out_row = [chrom, coord, ref_base, ref_counts, var_base, var_counts]
        
        writer.writerow(out_row)

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()

    parser.add_argument('vcf_file',
                        help='''Path to VCF file containing calls from the normal genome. This will be used to identify
                        SNPs.''')
    
    parser.add_argument('out_file',
                        help='''Path where output file will be written. The output format will be appropriate for input
                        into APOLLOH.''')
    
    parser.add_argument('--normal_column', type=int, default=0,
                        help='''0 based index of the column containing the normal sample. Choices are 0 or 1. The tumour
                        column will be assumed to be the other index.''')
    
    parser.add_argument('--min_genotype_quality', type=int, default=99)
    
    args = parser.parse_args()
    
    main(args)
