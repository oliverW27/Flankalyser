#!/usr/bin/env python3
import re
import argparse
from collections import defaultdict

def parse_tsv(tsv_file):
    """Parse lexicmap TSV to get qlen per query element."""
    qlen_lookup = {}
    with open(tsv_file, 'r') as f:
        header = f.readline().strip().split('\t')
        query_idx = header.index('query')
        qlen_idx = header.index('qlen')
        
        for line in f:
            cols = line.strip().split('\t')
            query = cols[query_idx]
            qlen = int(cols[qlen_idx])
            qlen_lookup[query] = qlen
    
    return qlen_lookup

def parse_fasta_header(header):
    """Extract all required fields from header."""
    biosample_match = re.search(r'sgenome=(\S+)', header)
    biosample = biosample_match.group(1) if biosample_match else 'unknown'
    
    query_match = re.search(r'query=(\S+)', header)
    query = query_match.group(1) if query_match else 'unknown'
    
    coords_match = re.search(r':(\d+)-(\d+):[+-]', header)
    seq_start = int(coords_match.group(1)) if coords_match else None
    seq_end = int(coords_match.group(2)) if coords_match else None
    
    sstart_match = re.search(r'sstart=(\d+)', header)
    send_match = re.search(r'send=(\d+)', header)
    sstart = int(sstart_match.group(1)) if sstart_match else None
    send = int(send_match.group(1)) if send_match else None
    
    qstart_match = re.search(r'qstart=(\d+)', header)
    qend_match = re.search(r'qend=(\d+)', header)
    qstart = int(qstart_match.group(1)) if qstart_match else None
    qend = int(qend_match.group(1)) if qend_match else None
    
    sstr_match = re.search(r'sstr=([+-])', header)
    sstr = sstr_match.group(1) if sstr_match else '+'
    
    return biosample, query, seq_start, seq_end, sstart, send, qstart, qend, sstr

def extract_flanks(sequence, seq_start, seq_end, sstart, send, qstart, qend, sstr, element_len, flank_size):
    """Extract left and right flanks based on hit coordinates and strand."""
    
    left_flank = None
    right_flank = None
    
    if sstr == '+':
        if qstart == 1:
            left_flank_len = sstart - seq_start
            if left_flank_len > 0:
                left_flank = sequence[:left_flank_len]
                if len(left_flank) < flank_size:
                    left_flank = 'N' * (flank_size - len(left_flank)) + left_flank
        
        if qend == element_len:
            right_flank_len = seq_end - send
            if right_flank_len > 0:
                right_flank = sequence[-right_flank_len:]
                if len(right_flank) < flank_size:
                    right_flank = right_flank + 'N' * (flank_size - len(right_flank))
    
    else:  # Minus strand
        if qstart == 1:
            left_flank_len = seq_end - send
            if left_flank_len > 0:
                left_flank = sequence[:left_flank_len]
                if len(left_flank) < flank_size:
                    left_flank = 'N' * (flank_size - len(left_flank)) + left_flank
        
        if qend == element_len:
            right_flank_len = sstart - seq_start
            if right_flank_len > 0:
                right_flank = sequence[-right_flank_len:]
                if len(right_flank) < flank_size:
                    right_flank = right_flank + 'N' * (flank_size - len(right_flank))
    
    return left_flank, right_flank

def main():
    parser = argparse.ArgumentParser(description="Extract flanking sequences from multifasta")
    parser.add_argument("--input", required=True, help="Input multifasta file")
    parser.add_argument("--tsv", required=True, help="Lexicmap TSV file with qlen info")
    parser.add_argument("--output_dir", required=True, help="Output directory for flank files")
    parser.add_argument("--flank_len", required=True, type=int, help="Flank length to extract")
    args = parser.parse_args()
    
    qlen_lookup = parse_tsv(args.tsv)
    query_data = defaultdict(lambda: defaultdict(lambda: {'left': [], 'right': []}))
    
    with open(args.input, 'r') as f:
        current_header = None
        current_seq = ""
        
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header and current_seq:
                    biosample, query, seq_start, seq_end, sstart, send, qstart, qend, sstr = parse_fasta_header(current_header)
                    element_len = qlen_lookup.get(query)
                    
                    if None not in (seq_start, seq_end, sstart, send, qstart, qend, element_len):
                        left, right = extract_flanks(
                            current_seq, seq_start, seq_end, sstart, send,
                            qstart, qend, sstr, element_len, args.flank_len
                        )
                        
                        if left is not None or right is not None:
                            query_data[query][biosample]['left'].append(left if left else 'N' * args.flank_len)
                            query_data[query][biosample]['right'].append(right if right else 'N' * args.flank_len)
                
                current_header = line[1:]
                current_seq = ""
            else:
                current_seq += line
        
        if current_header and current_seq:
            biosample, query, seq_start, seq_end, sstart, send, qstart, qend, sstr = parse_fasta_header(current_header)
            element_len = qlen_lookup.get(query)
            
            if None not in (seq_start, seq_end, sstart, send, qstart, qend, element_len):
                left, right = extract_flanks(
                    current_seq, seq_start, seq_end, sstart, send,
                    qstart, qend, sstr, element_len, args.flank_len
                )
                
                if left is not None or right is not None:
                    query_data[query][biosample]['left'].append(left if left else 'N' * args.flank_len)
                    query_data[query][biosample]['right'].append(right if right else 'N' * args.flank_len)
    
    for query, biosample_data in query_data.items():
        output_file = f"{args.output_dir}/flanks_{query}.txt"
        with open(output_file, 'w') as f:
            for biosample, data in biosample_data.items():
                left_str = '|'.join(data['left'])
                right_str = '|'.join(data['right'])
                f.write(f"{biosample}\t{left_str}\t{right_str}\n")
        print(f"Created: {output_file}")

if __name__ == "__main__":
    main()
