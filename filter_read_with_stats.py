#%%
from collections import Counter
import gzip, subprocess
import numpy as np
from zlib_ng import gzip_ng

#%%
def get_filesize_from_gzip(path_file):
    import subprocess
    
    cmd = f"gzip -dc {path_file}|wc -c"
    output = subprocess.check_output(cmd, shell = True)
    return int(output)

def estimate_uncompressed_gz_size(filename):
    from io import SEEK_END
    import io
    import os
    import struct
    import zlib
    import gzip
    # From the input file, get some data:
    # - the 32 LSB from the gzip stream
    # - 1MB sample of compressed data
    # - compressed file size
    with open(filename, "rb") as gz_in:
        gz_in.seek(-4, SEEK_END)
        lsb = struct.unpack('I', gz_in.read(4))[0]
        file_size = os.fstat(gz_in.fileno()).st_size
        
        
    read_bytes = 10000000
    
    with gzip.open(filename, "rb") as frb:
        sample = frb.read(read_bytes)
    
    compressed_sample = io.BytesIO()
    with gzip.GzipFile(fileobj = compressed_sample, mode = "wb") as fwb:
        fwb.write(sample)
    
    compressed_len = len(compressed_sample.getvalue())
    decompressed_len = len(sample)

    estimate = int(file_size * decompressed_len / compressed_len)

    # 32 LSB to zero
    mask = ~0xFFFFFFFF

    # Kill the 32 LSB to be substituted by the data read from the file
    adjusted_estimate = (estimate & mask) | lsb

    return adjusted_estimate

def get_kmer_count(sequence, nmer = 5):
    lenseq = len(sequence)
    
    list_kmer = list(map(lambda ind: sequence[ind:ind+nmer], range(lenseq-nmer+1)))
    return Counter(list_kmer)

def get_kmer_count_of_sequence_with_chop(sequence, chop_length, nmer):
    sequence_front = sequence[:chop_length]
    sequence_end = sequence[-chop_length:]
    sequence_intermediate = sequence[chop_length:-chop_length]
    
    count_kmer_front        = get_kmer_count(sequence_front, nmer)
    count_kmer_end          = get_kmer_count(sequence_end, nmer)
    count_kmer_intermediate = get_kmer_count(sequence_intermediate, nmer)
    
    return (count_kmer_front, count_kmer_end, count_kmer_intermediate)
    
class FASTA_READER:
    def __init__(self, path_fasta, is_gzip):
        self.is_gzip = is_gzip
        if is_gzip:
            self.file_reader = gzip_ng.open(path_fasta, "rb")
        else:
            self.file_reader = open(path_fasta, 'r')
    
    def readline(self):
        if self.is_gzip:
            return self.file_reader.readline().decode()
        else:
            return self.file_reader.readline()
    
    def seek(self, val):
        self.file_reader.seek(val)
    
    def tell(self):
        return self.file_reader.tell()

    def close(self):
        self.file_reader.close()
    
    def __del__(self):
        self.close()
        
class FASTA_WRITER:
    def __init__(self, path_fasta, is_gzip):
        self.is_gzip = is_gzip
        if is_gzip:
            self.file_reader = gzip_ng.open(path_fasta, "wb")
        else:
            self.file_reader = open(path_fasta, 'w')
    
    def write(self, value:str):
        if self.is_gzip:
            return self.file_reader.write(value.encode())
        else:
            return self.file_reader.write(value)
        
    def writeline(self, value:str):
        value_end = f"{value}\n"
        self.write(value_end)

    def close(self):
        self.file_reader.close()
    
    def __del__(self):
        self.close()
        
def calculate_mean_error_rate_as_phred_scale(arr_phredscaled_quality):
    arr_error_rate = 10 ** (-arr_phredscaled_quality/10)
    mean_error_rate = np.mean(arr_error_rate)
    mean_error_rate_phredscaled = -10 * np.log10(mean_error_rate)
    
    return mean_error_rate_phredscaled

def filter_stat_from_file_range(path_fasta, path_save_tmp, is_gzip, min_length, min_qv, min_q30_perc, readid_prefix = None, readid_suffix = None):
    import numpy as np
    
    filereader = FASTA_READER(path_fasta, is_gzip)
    filewriter = FASTA_WRITER(path_save_tmp, is_gzip)
    
    set_atgc = set("ATGCN")
    
    n_passed = 0
    n_filtered_by_length = 0
    n_filtered_by_qv = 0
    n_filtered_by_q30 = 0
    
    while 1:
        namevalue = filereader.readline()
        if namevalue == '':
            break
        if namevalue.startswith('@'):
            name = namevalue.strip("@\n")
            if readid_prefix != None:
                if not name.startswith(readid_prefix):
                    continue
            if readid_suffix != None:
                if not name.endswith(readid_suffix):
                    continue
        else:
            continue
        sequencevalue = filereader.readline()
        tmp = filereader.readline()
        qualityvalue = filereader.readline()
        
        sequence = sequencevalue.strip('\n').upper()
        quality = qualityvalue.strip('\n')
        
        assert set(sequence) <= set_atgc, name
        
        list_quality_int = list(map(lambda val: ord(val)-33, quality))
        
        qv = calculate_mean_error_rate_as_phred_scale(np.array(list_quality_int))
        len_read = len(quality)
        
        write_line = True
        if min_length:
            if min_length > len_read:
                write_line = False
                n_filtered_by_length += 1
        if min_qv:
            if min_qv > qv:
                if write_line:
                    n_filtered_by_qv += 1
                write_line = False
        if min_q30_perc:
            n_q30_base = sum(list(map(lambda val: val >= 30, list_quality_int)))
            q30_perc = n_q30_base / len_read * 100
            if min_q30_perc > q30_perc:
                if write_line:
                    n_filtered_by_q30 += 1
                write_line = False
        if write_line:
            n_passed += 1
            filewriter.write(namevalue)
            filewriter.write(sequencevalue)
            filewriter.write(tmp)
            filewriter.write(qualityvalue)
        
        n_processed = n_passed+n_filtered_by_length+n_filtered_by_qv+n_filtered_by_q30
        if (n_processed) % 10000 == 0:
            print(n_processed)
        
    filereader.close()
    filewriter.close()
    return (n_passed, n_filtered_by_length, n_filtered_by_qv, n_filtered_by_q30)

def filter_with_stat_from_file_with_ncores(path_fasta, path_save, min_length = None, min_qv = None, min_q30_perc = None, is_gzip = True, readid_prefix = None, readid_suffix = None):
    import os
    import numpy as np    
    
    os.makedirs(os.path.dirname(path_save), exist_ok=True)
    n_passed, n_filtered_by_length, n_filtered_by_qv, n_filtered_by_q30 = filter_stat_from_file_range(
        path_fasta, path_save, is_gzip, min_length, min_qv, min_q30_perc, readid_prefix, readid_suffix
    )
    
    with open(f"{path_save}.stat", 'w') as fw:
        fw.write(f"#Total : {n_passed+n_filtered_by_length+n_filtered_by_qv}\n")
        fw.write(f"#Passed : {n_passed}\n")
        fw.write(f"#Filtered_by_Length : {n_filtered_by_length}\n")
        fw.write(f"#Filtered_by_QV : {n_filtered_by_qv}\n")
        fw.write(f"#Filtered_by_Q30_Rate : {n_filtered_by_q30}\n")

    
#%%
if __name__ == "__main__":
    list_seqname = [
        "20240517_1158_P2S-01506-B_PAW01562_215af4e8",
        "20240629_1736_P2S-01506-A_PAW01687_bf35a486",
        "20240629_1736_P2S-01506-B_PAU99972_b4ed6da3",
        "20240811_1107_P2S-01506-A_PAW77551_0cd923fa",
        "20240816_0603_P2S-01501-B_PAW77472_d96fb6c9",
        "20240817_0645_P2S-01501-B_PAW77472_c8fb51b0",
        "20240818_1702_P2S-01506-B_PAW88380_5c7050f8",
        "20240920_1805_P2S-01506-A_PAY14781_7526b2f4",
        "20240923_1616_P2S-01506-A_PAY08433_74256187",
        "20240923_1641_P2S-01506-B_PAY03655_fb6f2a12",
        "20240924_2038_P2S-01506-A_PAY08433_7a23907e",
        "20240924_2038_P2S-01506-B_PAY03655_f02a3d27",
        "20240926_1205_P2S-01506-A_PAY09654_90f482e1",
        "20240926_1205_P2S-01506-B_PAY09817_4744eb86",
        "20240930_1539_P2S-01506-B_PAY10807_076f0b43",
        "20241002_1732_P2S-01506-A_PAY10807_3ab61abe",
        "20241011_1754_P2S-01506-A_PAW99288_c001806f",
        "20241012_1341_P2S-01506-A_PAW99288_c255cb62",
        "20241012_1341_P2S-01506-B_PAY01199_fef5a120",
        "20241021_1517_P2S-01506-A_PAW99288_e821305c",
        "20241021_1517_P2S-01506-B_PAY01199_b712d259",
        "20241022_1541_P2S-01506-A_PAW99013_6c642a53",
        "20241022_1541_P2S-01506-B_PAY07547_1caa91a9",
        "20241031_1212_P2S-01506-A_PAY83904_baf56bf4",
        "20241031_1220_P2S-01506-B_PAY74072_775e6f78",
        "20241031_1440_P2S-01503-A_PAY74208_5cb0c269",
        "20241031_1440_P2S-01503-B_PAY74099_bed3405c",
        "20241101_1326_P2S-01506-A_PAY83904_5bb2fcbf",
        "20241101_1326_P2S-01506-B_PAY74072_95afdfdf"
    ]
    from joblib import Parallel, delayed
    
    with Parallel(n_jobs=len(list_seqname)) as parallel:
        parallel(delayed(filter_with_stat_from_file_with_ncores)(
            f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/ONT/Edge_Trimmed/{seqname}.Adapt_Filt.Edge_Trim_Front_1500_Tail_1500.fastq",
            f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/ONT/UL_Filtered/{seqname}.Adapt_Filt.Edge_Trim_Front_1500_Tail_1500.Filter_L50kb_QV15_Q30Perc80.fastq",
            min_length = 50000,
            min_qv = 15,
            min_q30_perc = 80,
            is_gzip = False,
            readid_prefix = None,
            readid_suffix = "40GB"
        ) for seqname in list_seqname)