#%%
import gzip,re,json,os,sys
from collections import Counter

from zlib_ng import gzip_ng


def run(path_fastq, path_save, path_fastqids):
    list_fastqids_subseq = read_readids_subseq(path_fastqids)
    
    fastq_reader = get_fastq_reader(path_fastq)
    fastq_writer = get_fastq_writer(path_save)
    print(f"Start reading {path_fastq}...", flush = True)
    extract_read_by_fastq(fastq_reader, fastq_writer, list_fastqids_subseq)
    print(f"Writing {path_save} Finished", flush = True)
    close_buffer(fastq_reader)
    close_buffer(fastq_writer)
    
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
    
def extract_read_by_fastq(reader, writer, headers):
    set_headers = set(headers)
    
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
            write_read_set(writer, read_set)
            set_headers.remove(seqname)
            
            if len(set_headers) == 0:
                break
    print(f"{len(set_headers)} headers not found.", flush = True, file = sys.stdout)

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
    
#%%
if __name__ == "__main__":
    import argparse
    
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument("--input")
    argument_parser.add_argument("--out")
    argument_parser.add_argument("--file")
    
    args = argument_parser.parse_args()
    
    run(args.input, args.out, args.file)