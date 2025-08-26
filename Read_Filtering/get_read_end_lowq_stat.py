#%%
from collections import Counter
import gzip
import numpy as np

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
            self.file_reader = gzip.open(path_fasta, "rb")
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
        
def calculate_mean_error_rate_as_phred_scale(arr_phredscaled_quality):
    arr_error_rate = 10 ** (-arr_phredscaled_quality/10)
    mean_error_rate = np.mean(arr_error_rate)
    mean_error_rate_phredscaled = -10 * np.log10(mean_error_rate)
    
    return mean_error_rate_phredscaled

def get_stat_from_file_range(path_fasta, start, end, is_gzip, kmer_chop_length, n_kmer, index_low_quality, readid_prefix = None, readid_suffix = None):
    import numpy as np
    
    dict_stat = {
        "mean_quality_by_read":list(),
        "qv_by_read":list(),
        "std_quality_by_read":list(),
        "length_by_read":list(),
        "n_N_base":list(),
        "n_Q30_base":list(),
        "ind_lowQ_base":list(),
        "kmer_front":Counter(),
        "kmer_end":Counter(),
        "kmer_intermediate":Counter(),
        "readid":list()
    }
    
    filereader = FASTA_READER(path_fasta, is_gzip)
    filereader.seek(start)
    set_atgc = set("ATGCN")
    while 1:
        name = filereader.readline()
        if name == '':
            break
        if name.startswith('@'):
            name = name.strip("@\n")
            if readid_prefix != None:
                if not name.startswith(readid_prefix):
                    continue
            if readid_suffix != None:
                if not name.endswith(readid_suffix):
                    continue
        else:
            continue
        sequence = filereader.readline().strip('\n').upper()
        _ = filereader.readline()
        quality = filereader.readline().strip('\n')
        
        assert set(sequence) <= set_atgc, name
        
        list_quality_int = list(map(lambda val: ord(val)-33, quality))
        n_q30_base = sum(list(map(lambda val: val >= 30, list_quality_int)))
        
        dict_stat["mean_quality_by_read"].append(np.mean(list_quality_int))
        dict_stat["qv_by_read"].append(calculate_mean_error_rate_as_phred_scale(np.array(list_quality_int)))
        dict_stat["std_quality_by_read"].append(np.std(list_quality_int))
        dict_stat["length_by_read"].append(len(quality))
        dict_stat["n_N_base"].append(Counter(sequence)['N'])
        dict_stat["n_Q30_base"].append(n_q30_base)
        dict_stat["readid"].append(name)
        
        if index_low_quality != None:
            dict_stat["ind_lowQ_base"].append(
                tuple(filter(lambda ind: list_quality_int[ind] < index_low_quality, range(len(list_quality_int))))
            )
        
        if n_kmer != None:
            count_kmer_front, count_kmer_end, count_kmer_intermediate = get_kmer_count_of_sequence_with_chop(sequence, kmer_chop_length, n_kmer)
            for name, cnt in count_kmer_front.items():
                dict_stat["kmer_front"][name] += cnt
            for name, cnt in count_kmer_end.items():
                dict_stat["kmer_end"][name] += cnt
            for name, cnt in count_kmer_intermediate.items():
                dict_stat["kmer_intermediate"][name] += cnt
        
        current_cursor = filereader.tell()
        if end > 0 and current_cursor > end:
            break
    filereader.close()
    return dict_stat

def get_stat_from_file_with_ncores(path_fasta, n_cpu, kmer_chop_length, n_kmer, index_low_quality = None, filesize = None, is_gzip = True, readid_prefix = None, readid_suffix = None):
    import os
    import numpy as np
    from joblib import Parallel, delayed
    
    if filesize == None:
        if is_gzip:
            filesize = estimate_uncompressed_gz_size(path_fasta)
        else:
            filesize = os.stat(path_fasta).st_size
    
    list_ranges = list(map(int, np.linspace(0, filesize, n_cpu+1)))
    list_ranges[-1] = -1
    
    list_starts = list_ranges[:-1]
    list_ends = list_ranges[1:]
    
    with Parallel(n_jobs=n_cpu) as parallel:
        list_stats = parallel(delayed(get_stat_from_file_range)(
            path_fasta, start, end, is_gzip, kmer_chop_length, n_kmer, index_low_quality, readid_prefix, readid_suffix
        )for start, end in zip(list_starts, list_ends))
    
    dict_stat = {
        "mean_quality_by_read":list(),
        "qv_by_read":list(),
        "std_quality_by_read":list(),
        "length_by_read":list(),
        "n_N_base":list(),
        "n_Q30_base":list(),
        "ind_lowQ_base":list(),
        "kmer_front":Counter(),
        "kmer_end":Counter(),
        "kmer_intermediate":Counter(),
        "readid":list()
    }
    
    for stat_tmp in list_stats:
        for key in dict_stat.keys():
            if isinstance(dict_stat[key], list):
                dict_stat[key].extend(stat_tmp[key])
            if isinstance(dict_stat[key], Counter):
                for name, cnt in stat_tmp[key].items():
                    dict_stat[key][name] += cnt
    return dict_stat

def get_lowquality_base_count_from_end(dict_stat, check_length):
    import numpy as np
    cnt_front = np.zeros(check_length)
    cnt_end = np.zeros(check_length)
    
    checkarange = np.arange(check_length)
    
    for length, list_lowqual_ind in zip(dict_stat["length_by_read"], dict_stat["ind_lowQ_base"]):
        arr_lowqual_ind_front = np.array(list_lowqual_ind)
        arr_lowqual_ind_end = (length - 1) - arr_lowqual_ind_front
        
        cnt_front += np.isin(checkarange, arr_lowqual_ind_front[:check_length])
        cnt_end += np.isin(checkarange, arr_lowqual_ind_end[:check_length])
        
    return cnt_front, cnt_end

#%%
import pickle

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

# %%
################################
dict_seqname_to_lowqcnt_untrimmed = dict()
dict_seqname_to_lowqcnt_trimmed = dict()
#%%
import pickle
import numpy as np

for seqname in list_seqname:
    list_name_val = seqname.split('_')
    seqname_repr = f"{list_name_val[0]}:{list_name_val[-2]}"
    
    with open(f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/ONT/Adapter_Filtered/Stat_pickle/Stat_{seqname}.Adapt_Filt.pickle", "rb") as frb:
        dict_seqname_to_stat_untrimmed = pickle.load(frb)
        
    with open(f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/ONT/Edge_Trimmed/Stat_pickle/Stat_{seqname}.Adapt_Filt.Edge_Trim_Front_1500_Tail_1500.pickle", "rb") as frb:
        dict_seqname_to_stat_trimmed = pickle.load(frb)
        
    cnt_front_unt, cnt_end_unt = get_lowquality_base_count_from_end(dict_seqname_to_stat_untrimmed, 5000)
    cnt_front_t, cnt_end_t = get_lowquality_base_count_from_end(dict_seqname_to_stat_trimmed, 5000)
    
    n_read_unt = len(dict_seqname_to_stat_untrimmed["readid"])
    mean_perc_lowquality_unt = sum(list(map(len, dict_seqname_to_stat_untrimmed["ind_lowQ_base"]))) / sum(dict_seqname_to_stat_untrimmed["length_by_read"]) * 100
    
    n_read_t = len(dict_seqname_to_stat_trimmed["readid"])
    mean_perc_lowquality_t = sum(list(map(len, dict_seqname_to_stat_trimmed["ind_lowQ_base"]))) / sum(dict_seqname_to_stat_trimmed["length_by_read"]) * 100
    
    dict_seqname_to_lowqcnt_untrimmed[seqname_repr] = {
        "cnt_front" : cnt_front_unt,
        "cnt_end" : cnt_end_unt,
        "n_read" : n_read_unt,
        "mean_perc_lowquality" : mean_perc_lowquality_unt
    }
    dict_seqname_to_lowqcnt_trimmed[seqname_repr] = {
        "cnt_front" : cnt_front_t,
        "cnt_end" : cnt_end_t,
        "n_read" : n_read_t,
        "mean_perc_lowquality" : mean_perc_lowquality_t
    }
    
    del dict_seqname_to_stat_trimmed, dict_seqname_to_stat_untrimmed
    print(seqname)

with open(f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/ONT/Adapter_Filtered/Stat_pickle/Edge_Lowq_Stat.Adapt_Filt.pickle", "wb") as fwb:
    pickle.dump(dict_seqname_to_lowqcnt_untrimmed, fwb)
    
with open(f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/ONT/Edge_Trimmed/Stat_pickle/Edge_Lowq_Stat.Adapt_Filt.Edge_Trim_Front_1500_Tail_1500.pickle", "wb") as fwb:
    pickle.dump(dict_seqname_to_lowqcnt_trimmed, fwb)