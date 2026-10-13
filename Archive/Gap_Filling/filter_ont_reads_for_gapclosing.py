#%%
import gzip,re,json,os,sys
from collections import Counter

from zlib_ng import gzip_ng
import pandas as pd
import numpy as np


def run(path_fastq, path_save, path_mapped_paf, qry_align_ratio, identity_ratio, aligned_region):
    table_paf = read_paf(path_mapped_paf)
    list_perfect_matched_reads = get_perfect_matched_reads(table_paf, qry_align_ratio, identity_ratio, aligned_region)
    save_fastq_without_perfectmatch(path_fastq, path_save, list_perfect_matched_reads)

def read_paf(path_paf):
    table_paf = pd.read_csv(path_paf, sep = '\t', names = list(range(1, 13)), usecols = list(range(1, 13)))
    table_paf = table_paf.rename(
        columns = {
            1:"qryname",
            2:"qrylen",
            3:"qry_alignstart",
            4:"qry_alignend",
            5:"strand",
            6:"refname",
            7:"reflen",
            8:"ref_alignstart",
            9:"ref_alignend",
            10:"match",
            11:"align",
            12:"quality"
        }
    )
    return table_paf

def get_perfect_matched_reads(table_paf, qry_align_ratio, identity_ratio, aligned_region):
    table_paf_high_qual = table_paf[table_paf["quality"] >= 60]
    table_paf_high_qual["qry_alignlen"] = table_paf_high_qual["qry_alignend"] - table_paf_high_qual["qry_alignstart"]
    table_paf_high_qual["qry_align_ratio"] = table_paf_high_qual["qry_alignlen"] / table_paf_high_qual["qrylen"]
    table_paf_high_qual["match_to_align_ratio"] = table_paf_high_qual["match"] / table_paf_high_qual["align"]

    table_paf_high_qual["align_inrange_start"] = table_paf_high_qual["ref_alignstart"] > aligned_region
    table_paf_high_qual["align_inrange_end"] = table_paf_high_qual["ref_alignend"] < (table_paf_high_qual["reflen"] - aligned_region)

    table_paf_high_qual_ratio = table_paf_high_qual[table_paf_high_qual["qry_align_ratio"] > qry_align_ratio]

    table_paf_high_qual_ratio_match = table_paf_high_qual_ratio[table_paf_high_qual_ratio["match_to_align_ratio"] > identity_ratio]

    table_paf_high_qual_ratio_match_inregion = table_paf_high_qual_ratio_match[np.logical_and(table_paf_high_qual_ratio_match["align_inrange_start"], table_paf_high_qual_ratio_match["align_inrange_end"])]

    return table_paf_high_qual_ratio_match_inregion["qryname"].to_list()

def get_fastq_reader(path_file):
    if not os.path.exists(path_file):
        return iter([])
    if path_file.endswith("gz"):
        reader = gzip_ng.open(path_file, 'rb')
    else:
        reader = open(path_file, 'r')
    return reader

def get_fastq_writer(path_file):
    if path_file.endswith("gz"):
        writer = gzip_ng.open(path_file, 'wb')
    else:
        writer = open(path_file, 'w')
    return writer

def read_readids_subseq(path_file):
    with open(path_file, 'r') as fr:
        list_lines = fr.readlines()
    
    list_readids = list(map(lambda val: val.strip('\n'), list_lines))
    return list_readids

def close_buffer(file_buffer):
    if hasattr(file_buffer, "closed") and not file_buffer.closed:
        file_buffer.close()

def save_fastq_without_perfectmatch(path_fastq, path_save, list_perfect_matched_reads):
    fastq_reader = get_fastq_reader(path_fastq)
    fastq_writer = get_fastq_writer(path_save)
    extract_fastq_except_reads(fastq_reader, fastq_writer, list_perfect_matched_reads)
    
    close_buffer(fastq_reader)
    close_buffer(fastq_writer)

def extract_fastq_except_reads(reader, writer, headers):
    set_headers = set(headers)
    list_headers_found = list()

    while True:
        read_set = get_single_read_set(reader)
        is_set_valid, is_set_end = check_read_set_valid(read_set)
        if is_set_end: break
        if not is_set_valid:
            print(f"{reader.name} File has unrecognizable read. Pass this read instead.\n\t{read_set[0]}\n\t{read_set[1]}\n\t{read_set[2]}\n\t{read_set[3]}\n", flush = True, file = sys.stderr)
            continue
        
        assert read_set[0].startswith('@'), read_set
        seqname = read_set[0].strip('\n')[1:].split()[0]

        if seqname in set_headers:
            list_headers_found.append(seqname)
        else:
            write_read_set(writer, read_set)
            
    print(f"{len(list_headers_found)} headers skipped.", flush = True, file = sys.stdout)
    print(f"{len(set_headers - set(list_headers_found))} headers out of {len(set_headers)} not found.", flush = True, file = sys.stdout)

def get_single_read_set(reader):
    seqid = read_line_from_reader(reader, None)
    seq = read_line_from_reader(reader, None)
    optional = read_line_from_reader(reader, None)
    quality = read_line_from_reader(reader, None)
    return (seqid, seq, optional, quality)

def check_read_set_valid(read_set):
    list_is_valid = list(map(lambda x : x != None, read_set))
    is_set_valid = (sum(list_is_valid) == 4)
    is_set_end = (sum(list_is_valid) == 0)
    return is_set_valid, is_set_end

def read_line_from_reader(reader, replace):
    line = next(reader, replace)
    if line != None and hasattr(line, "decode"):
        line = line.decode()
    return line   

def write_read_set(writer, read_set):
    for line in read_set:
        write_line_by_writer(writer, line)

def write_line_by_writer(writer, line):
    line = line.strip('\n') + '\n'
    if hasattr(writer, "compress"):
        line = line.encode()
    writer.write(line)
# %%
if __name__ == '__main__':
    import argparse
    
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument("--fastq_in")
    argument_parser.add_argument("--fastq_out")
    argument_parser.add_argument("--paf")
    argument_parser.add_argument("--align_ratio", type = float, default = 0.95)
    argument_parser.add_argument("--identity_ratio", type = float, default = 0.8)
    argument_parser.add_argument("--aligned_region", type = int, default = 1000000)
    
    args = argument_parser.parse_args()

    run(args.fastq_in, args.fastq_out, args.paf, args.align_ratio, args.identity_ratio, args.aligned_region)
# %%
