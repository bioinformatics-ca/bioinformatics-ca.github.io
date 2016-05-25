import csv

def main(args):
    states = load_states(args.state_file)
    
    segments = load_segments(args.cnv_file, states)
    
    write_copy_number_file(segments, args.out_file)

def load_states(file_name):
    states = {}
    
    reader = csv.DictReader(open(file_name), delimiter='\t')
    
    for i, row in enumerate(reader):
        state = i + 1
        
        states[state] = {}
        
        states[state]['major_cn'] = int(row['TumourBalleles3'])
        
        states[state]['minor_cn'] = int(row['TumourBalleles2'])
    
    return states

def load_segments(file_name, states):
    segments = []
    
    reader = csv.DictReader(open(file_name), delimiter='\t')
    
    for row in reader:
        out_row = {}
        
        chrom = row['Chromosome']
        
        if chrom == '23':
            chrom = 'X'
        elif chrom == '24':
            chrom = 'Y'
        
        out_row['chrom'] = chrom 
        
        out_row['beg'] = row['Start Position (bp)']
        
        out_row['end'] = row['End Position (bp)']
        
        segment_state = int(row['Tumour State']) 
        
        out_row['major_cn'] = states[segment_state]['major_cn']
        
        out_row['minor_cn'] = states[segment_state]['minor_cn']
        
        out_row['total_cn'] = int(row['Copy Number'])
        
        segments.append(out_row)
    
    return segments

def write_copy_number_file(segments, out_file):
    fields = ['chrom', 'beg', 'end', 'major_cn', 'minor_cn', 'total_cn']
    
    writer = csv.DictWriter(open(out_file, 'w'), fields, delimiter='\t')
    
    writer.writeheader()
    
    writer.writerows(segments)

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('cnv_file', help='''Path to .cnv file produced by OncoSNP.''')
    
    parser.add_argument('out_file', help='''Path where output file will be written in tsv format.''')
    
    parser.add_argument('state_file', help='''OncoSNP state file.''')
    
    args = parser.parse_args()
    
    main(args)
