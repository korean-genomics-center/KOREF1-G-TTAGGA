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
# %%
n_kmer = None
kmer_chop_length = None
index_lowquality = None
readid_suffix = "40GB"
#%%
path_fasta = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/ONT/Merged/ONT_UL.Adapt_Filt.Edge_Trim_Front_1500_Tail_1500.Filter_L50kb_QV20_Q30Perc80.merged.fastq"

dict_stat = get_stat_from_file_with_ncores(path_fasta, n_cpu=100, kmer_chop_length=kmer_chop_length, n_kmer=n_kmer, index_low_quality = index_lowquality, filesize = None, is_gzip = False, readid_suffix = readid_suffix)

#%%
import pickle

with open("/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/ONT/Adapter_Filtered/Stat_pickle/Stat_20240517_1158_P2S-01506-B_PAW01562_215af4e8.Adapt_Filt.pickle", "rb") as frb:
    dict_stat = pickle.load(frb)

#%%
def make_kmer_count_dataframe(count_kmer, n_kmer = 5, contexts = "ATCG"):
    from itertools import product
    import pandas as pd
    import math
    
    n_mer_index = math.ceil(n_kmer/2)
    n_mer_column = n_kmer - n_mer_index
    
    list_context_index = list(map(lambda vals: ''.join(vals), product(contexts, repeat = n_mer_index)))
    list_context_column = list(map(lambda vals: ''.join(vals), product(contexts, repeat = n_mer_column)))
    
    table_kmer_count = pd.DataFrame(index = list_context_index, columns = list_context_column).fillna(0).astype(int)
    
    for context, cnt in count_kmer.items():
        context_ind = context[:n_mer_index]
        context_col = context[-n_mer_column:]
        
        table_kmer_count.loc[context_ind, context_col] = cnt
    
    return table_kmer_count
    
# %%
table_kmer_front = make_kmer_count_dataframe(
    dict_stat["kmer_front"], n_kmer
)
table_kmer_end = make_kmer_count_dataframe(
    dict_stat["kmer_end"], n_kmer
)
table_kmer_intermediate = make_kmer_count_dataframe(
    dict_stat["kmer_intermediate"], n_kmer
)
# %%
import seaborn as sns
from matplotlib import pyplot as plt
plt.rcParams["font.size"] = 10

plt.figure(figsize = (15, 13))
sns.heatmap(data = table_kmer_front, cmap = "flare_r")
plt.show()
# %%
plt.figure(figsize = (15, 13))
sns.heatmap(data = table_kmer_end, cmap = "flare_r")
plt.show()
#%%
plt.figure(figsize = (15, 13))
sns.heatmap(data = table_kmer_intermediate, cmap = "flare_r")
plt.show()
# %%
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

def draw_length_to_quality_scatter(dict_stat):
    import math
    
    plt.rcParams["font.size"] = 14
    fig, ax = plt.subplots(
        2, 2, gridspec_kw = {"width_ratios":(3, 1), "height_ratios":(1, 3)}, figsize = (7, 7)
    )
    plt.subplots_adjust(
        wspace = 0.15, hspace = 0.15
    )
    
    max_length = max(dict_stat["length_by_read"])
    length_showlim = math.ceil(max_length * 1.05 / 1000) * 1000
    
    ax[0, 0].hist(dict_stat["length_by_read"], bins = np.linspace(0, length_showlim, 51), color = "gray", zorder = 1)
    ax[0, 0].set_xlim(0, length_showlim)
    # ax[0, 0].set_xticks([])
    ax[0, 0].set_ylabel("# Reads")
    ax[0, 0].axvline(np.mean(dict_stat["length_by_read"]), linewidth = 2, linestyle = "--", color = "firebrick", zorder = 3)

    ax[0, 1].axis("off")

    ax[1, 0].scatter(
        dict_stat["length_by_read"],
        dict_stat["qv_by_read"],
        s = 5,
        color = "gray"
    )
    ax[1, 0].set_xlim(0, length_showlim)
    ax[1, 0].set_ylim(0, 95)
    ax[1, 0].set_xlabel("Read Length (bp)")
    ax[1, 0].set_ylabel("Quality Value (QV, Phred-scaled)")

    ax[1, 1].hist(dict_stat["qv_by_read"], bins = np.linspace(0, 95, 51), color = "gray", zorder = 1, orientation = "horizontal")
    ax[1, 1].set_ylim(0, 95)
    # ax[1, 1].set_yticks([])
    ax[1, 1].set_xlabel("# Reads")
    ax[1, 1].axhline(np.mean(dict_stat["qv_by_read"]), linewidth = 2, linestyle = "--", color = "firebrick", zorder = 3)

    plt.show()
#%%
def print_stat(dict_stat):
    import numpy as np
    
    print(len(dict_stat["readid"]))
    print(sum(dict_stat["length_by_read"]))
    print(sum(dict_stat["n_N_base"]))
    print()
    print(np.mean(dict_stat["length_by_read"]))
    print(min(dict_stat["length_by_read"]))
    print(np.quantile(dict_stat["length_by_read"], 0.01))
    print(np.quantile(dict_stat["length_by_read"], 0.99))
    print(max(dict_stat["length_by_read"]))
    print()
    print(sum(dict_stat["n_Q30_base"]) / sum(dict_stat["length_by_read"]) * 100)
    print(np.mean(dict_stat["qv_by_read"]))
    print(sum(np.array(dict_stat["mean_quality_by_read"]) * np.array(dict_stat["length_by_read"])) / sum(dict_stat["length_by_read"]))
    print(-10 * np.log10(sum((10 ** (-np.array(dict_stat["qv_by_read"]) / 10)) * np.array(dict_stat["length_by_read"])) / sum(dict_stat["length_by_read"])))

#%%
import pandas as pd
from matplotlib import pyplot as plt

def plot_quality_scatterplot(dict_stat):
    table_stat = pd.DataFrame(
        {
            key : dict_stat[key] for key in ['mean_quality_by_read', 'std_quality_by_read', 'qv_by_read', 'length_by_read', 'n_N_base', 'n_Q30_base', 'n_LowQ_base', 'readid']
        }    
    )
    table_stat["perc_Q30_base"] = table_stat["n_Q30_base"] / table_stat["length_by_read"] * 100
    table_stat["perc_lowQ_base"] = table_stat["n_LowQ_base"] / table_stat["length_by_read"] * 100
    
    dict_key_to_show = {
        "mean_quality_by_read" : "Mean Base Phred Quality Score",
        "qv_by_read" : "Quality Value per Read",
        # "std_quality_by_read" : "STD Base Quality",
        "perc_Q30_base" : "Q30 Base (%)",
        "perc_lowQ_base" : "<Q10 Base (%)",
        "length_by_read" : "Read Length (bp)"
    }
    
    fig, ax = plt.subplots(
        len(dict_key_to_show), len(dict_key_to_show), figsize = (4 * len(dict_key_to_show), 4 * len(dict_key_to_show))
    )
    plt.subplots_adjust(
        wspace = 0.25, hspace = 0.25
    )
    
    dict_item_to_lim = dict()
    for item in dict_key_to_show.keys():
        valmin = min(table_stat[item])
        valmax = max(table_stat[item])
        valwidth = valmax - valmin
        valadjust = valwidth / 20
        
        dict_item_to_lim[item] = (valmin - valadjust, valmax + valadjust)
    
    for ind_item1, item1 in enumerate(dict_key_to_show.keys()):
        for ind_item2, item2 in enumerate(dict_key_to_show.keys()):
            ax_here = ax[ind_item1, ind_item2]
            if ind_item2 == 0:
                ax_here.set_ylabel(dict_key_to_show[item1])
            if ind_item1 == len(dict_key_to_show)-1:
                ax_here.set_xlabel(dict_key_to_show[item2])
            
            if ind_item1 < ind_item2:
                ax_here.axis("off")
                continue
            if ind_item1 == ind_item2:
                lowerbound, upperbound = dict_item_to_lim[item1]
                ax_here.hist(
                    table_stat[item1],
                    color = "gray",
                    bins = np.linspace(lowerbound, upperbound, 67)
                )    
                ax_here.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
                ax_here.set_xlim(lowerbound, upperbound)
                # ax_here.set_ylabel("# Reads")
            else:
                xlim_lowerbound, xlim_upperbound = dict_item_to_lim[item2]
                ylim_lowerbound, ylim_upperbound = dict_item_to_lim[item1]
                
                ax_here.scatter(
                    table_stat[item2],
                    table_stat[item1],
                    s = 5,
                    color = "gray"
                )
                ax_here.set_xlim(xlim_lowerbound, xlim_upperbound)
                ax_here.set_ylim(ylim_lowerbound, ylim_upperbound)
                if item1 == "length_by_read":
                    ax_here.set_yticklabels([f"{ytickval/1000:.1f}k" for ytickval in ax_here.get_yticks()])
            # if item2 == "qv_by_read":
            #     ax_here.axvline(20, linewidth = 2, color = "firebrick", linestyle = '--')
            # elif item1 == "qv_by_read":
            #     ax_here.axhline(20, linewidth = 2, color = "firebrick", linestyle = '--')
    plt.show()
    
# %%
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
    
    

# %%
def draw_lowquality_perc_on_end(cnt_front, cnt_end, dict_stat, check_length = 2000):
    n_read = len(dict_stat["readid"])
    
    mean_perc_lowquality = sum(list(map(len, dict_stat["ind_lowQ_base"]))) / sum(dict_stat["length_by_read"]) * 100

    plt.rcParams["font.size"] = 14
    fig, ax = plt.subplots(1, 2, figsize = (12, 4))
    ax[0].plot(range(check_length), cnt_front / n_read * 100, color = "gray")
    ax[0].axhline(mean_perc_lowquality, linewidth = 1.5, linestyle = '--', color = "firebrick")
    ax[0].set_xlabel("Coordinate (bp)")
    ax[0].set_ylabel("Percentage Low Quality Base (%)")
    ax[0].set_title("Front")

    ax[1].plot(range(check_length), cnt_end / n_read * 100, color = "gray")
    ax[1].axhline(mean_perc_lowquality, linewidth = 1.5, linestyle = '--', color = "firebrick")
    ax[1].set_xlabel("Coordinate (bp)")
    ax[1].set_ylabel("Percentage Low Quality Base (%)")
    ax[1].set_title("End")
    plt.show()
# %%
######################################
import pickle
#%%
list_seqname = [
"20240510_1615_P2S-01506-A_PAU90966_ad72ac81",
"20240510_1615_P2S-01506-B_PAU91096_39bb4496",
"20240603_1341_P2S-01506-A_PAU69619_32394ac4",
"20240603_1341_P2S-01506-B_PAU69571_419b8c25",
"20240603_1553_P2S-01506-A_PAU69619_a90db93e",
"20240603_1553_P2S-01506-B_PAU69571_1464b73b",
"20240603_1925_P2S-01506-A_PAU69619_b2756a99",
"20240603_1925_P2S-01506-B_PAU69571_bd44489b",
"20240603_2030_P2S-01506-A_PAU69619_fad37a7f",
"20240603_2030_P2S-01506-B_PAU69571_54256320",
"20240604_1338_P2S-01506-A_PAU90972_437fd6b6",
"20240604_1338_P2S-01506-B_PAU91173_a1c4c20d",
"20240606_1426_P2S-01506-A_PAU91084_50a9eb35",
"20240606_1426_P2S-01506-B_PAU91079_79717851",
"20240608_1349_P2S-01506-A_PAU91135_57288ce9",
"20240608_1349_P2S-01506-B_PAU68508_19aa0ef8",
"20240611_1318_P2S-01506-A_PAU91094_93be457b"
]
#%%
dict_seqname_to_stat_duplex = dict()

for seqname in list_seqname:
    list_name_val = seqname.split('_')
    
    seqname_repr = f"{list_name_val[0]}:{list_name_val[-2]}:{list_name_val[-1]}"
    
    assert dict_seqname_to_stat_duplex.get(seqname_repr) == None
    with open(f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/ONT/Adapter_Filtered/Stat_pickle/Stat_{seqname}_duplex.pickle", "rb") as frb:
        dict_seqname_to_stat_duplex[seqname_repr] = pickle.load(frb)
        
dict_seqname_to_stat_simplex = dict()

for seqname in list_seqname:
    list_name_val = seqname.split('_')
    
    seqname_repr = f"{list_name_val[0]}:{list_name_val[-2]}:{list_name_val[-1]}"
    
    assert dict_seqname_to_stat_simplex.get(seqname_repr) == None
    with open(f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/ONT/Adapter_Filtered/Stat_pickle/Stat_{seqname}_simplex.pickle", "rb") as frb:
        dict_seqname_to_stat_simplex[seqname_repr] = pickle.load(frb)
        
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
#%%
dict_seqname_to_stat = dict()

for seqname in list_seqname:
    list_name_val = seqname.split('_')
    
    seqname_repr = f"{list_name_val[0]}:{list_name_val[-2]}"
    
    assert dict_seqname_to_stat.get(seqname_repr) == None
    with open(f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/ONT/Adapter_Filtered/Stat_pickle/Stat_{seqname}.Adapt_Filt.NoIndexing.pickle", "rb") as frb:
        dict_seqname_to_stat[seqname_repr] = pickle.load(frb)
#%%
dict_seqname_to_stat_unfilt = dict()

for seqname in list_seqname:
    list_name_val = seqname.split('_')
    
    seqname_repr = f"{list_name_val[0]}:{list_name_val[-2]}"
    
    assert dict_seqname_to_stat_unfilt.get(seqname_repr) == None
    with open(f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/ONT/Adapter_Filtered/Stat_pickle/Stat_{seqname}.pickle", "rb") as frb:
        dict_seqname_to_stat_unfilt[seqname_repr] = pickle.load(frb)
#%%
dict_seqname_to_stat_trimmed = dict()

for seqname in list_seqname:
    list_name_val = seqname.split('_')
    
    seqname_repr = f"{list_name_val[0]}:{list_name_val[-2]}"
    
    assert dict_seqname_to_stat_trimmed.get(seqname_repr) == None
    with open(f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/ONT/Edge_Trimmed/Stat_pickle/Stat_{seqname}.Adapt_Filt.Edge_Trim_Front_1500_Tail_1500.NoIndexing.pickle", "rb") as frb:
        dict_seqname_to_stat_trimmed[seqname_repr] = pickle.load(frb)
# %%
def mean_confidence_interval(data, confidence=0.95):
    from scipy import stats
    import numpy as np
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), stats.sem(a)
    h = se * stats.t.ppf((1 + confidence) / 2., n-1)
    return h


# %%
from matplotlib import pyplot as plt

plt.rcParams["font.size"] = 14
plt.figure(figsize = (10, 3))

xticks = list(dict_seqname_to_stat_simplex.keys())

plt.bar(list(map(lambda val: val - 0.15, range(len(xticks)))),
        list(map(lambda key: np.mean(dict_seqname_to_stat_simplex[key]["qv_by_read"]), xticks)),
        color = "gray",
        width = 0.3,
        label = "Simplex")
plt.bar(list(map(lambda val: val + 0.15, range(len(xticks)))),
        list(map(lambda key: np.mean(dict_seqname_to_stat_duplex[key]["qv_by_read"]), xticks)),
        color = "firebrick",
        width = 0.3, label = "Duplex")


plt.xticks(range(len(xticks)), xticks, rotation = -45, ha = "left")
plt.ylabel("Mean QV")
plt.legend(loc = "center left", bbox_to_anchor = (1.01, 0.5))
plt.show()
#%%
#%%
fig, ax = plt.subplots(4, 5, figsize = (12, 8))
axes = ax.flatten()
for ind, key in enumerate(dict_seqname_to_stat_duplex.keys()):
    
    perc_q30 = np.array(dict_seqname_to_stat_duplex[key]["n_Q30_base"]) / np.array(dict_seqname_to_stat_duplex[key]["length_by_read"]) * 100
    
    axes[ind].scatter(dict_seqname_to_stat_duplex[key]["qv_by_read"], perc_q30, color = "gray", s = 7)
    axes[ind].set_xlim([0, 50])
    axes[ind].set_ylim([0, 100])
    if ind < 15:
        pass
    else:
        # axes[ind].set_xticks([0, 15, 50])
        axes[ind].set_xlabel("QV per Read")
    if ind % 5 == 0:
        axes[ind].set_ylabel("Q30 Base (%)")
    # axes[ind].set_title(key, fontsize = 10)
    # axes[ind].ticklabel_format(axis = 'y', style = "sci", scilimits = (0, 0))
    # axes[ind].axvline(15, linewidth = 1.5, color = "firebrick", linestyle = "--")
# axes[-1].axis("off")

plt.tight_layout()
plt.show()
# %%
from matplotlib import pyplot as plt
import numpy as np

plt.rcParams["font.size"] = 14
plt.figure(figsize = (10, 3))

list_confidence_interval = list(map(lambda key: mean_confidence_interval(dict_seqname_to_stat[key]["qv_by_read"], 0.95), dict_seqname_to_stat.keys()))

plt.bar(dict_seqname_to_stat.keys(),
        list(map(lambda key: np.mean(dict_seqname_to_stat[key]["qv_by_read"]), dict_seqname_to_stat.keys())),
        yerr = list_confidence_interval,
        color = "gray")
plt.xticks(rotation = -45, ha = "left")
plt.ylabel("Mean QV")
plt.show()
# %%
plt.figure(figsize = (10, 3))
plt.bar(dict_seqname_to_stat.keys(),
        list(map(lambda key: len(dict_seqname_to_stat[key]["readid"]), dict_seqname_to_stat.keys())),
        color = "gray",
        label = "Total Reads",
        zorder = 0)

list_n_reads_over_length = list(map(lambda key:
                            len(list(filter(lambda length: length >= 100000, dict_seqname_to_stat[key]["length_by_read"])))
                            , dict_seqname_to_stat.keys()))
plt.bar(
    dict_seqname_to_stat.keys(),
    list_n_reads_over_length,
    color = "firebrick",
    label = "UL Reads (Length ≥ 100kb)",
    zorder = 3
)
for ind, key in enumerate(dict_seqname_to_stat.keys()):
    n_reads = len(dict_seqname_to_stat[key]["readid"])
    n_ul_reads = list_n_reads_over_length[ind]
    plt.annotate(f"{n_ul_reads/n_reads*100:.2f} %", (ind, n_reads+1e5), rotation = 90, va = "bottom", ha = "center")
plt.xticks(rotation = -45, ha = "left")
plt.ylabel("# Reads")
# plt.yscale("log")
plt.legend(loc = "center left", bbox_to_anchor = (1.01, 0.5))
plt.ylim(top = 6.5e6)
plt.show()
# %%
plt.figure(figsize = (10, 3))
plt.bar(dict_seqname_to_stat.keys(),
        list(map(lambda key: sum(dict_seqname_to_stat[key]["length_by_read"]), dict_seqname_to_stat.keys())),
        color = "gray",
        label = "Total Bases",
        zorder = 0)

list_ul_bases = list(map(lambda key: 
                            sum(list(filter(lambda length: length >= 100000, dict_seqname_to_stat[key]["length_by_read"])))
                            , dict_seqname_to_stat.keys()))

plt.bar(
    dict_seqname_to_stat.keys(),
    list_ul_bases,
    color = "firebrick",
    label = "UL Bases (Length ≥ 100kb)",
    zorder = 3
)
for ind, key in enumerate(dict_seqname_to_stat.keys()):
    bases = sum(dict_seqname_to_stat[key]["length_by_read"])
    n_ul_bases = list_ul_bases[ind]
    plt.annotate(f"{n_ul_bases/bases*100:.2f} %", (ind, bases+1e9), rotation = 90, va = "bottom", ha = "center")
plt.xticks(rotation = -45, ha = "left")
plt.ylabel("Bases (bp)")
# plt.yscale("log")
plt.legend(loc = "center left", bbox_to_anchor = (1.01, 0.5))
plt.ylim(top = 1.4e11)
plt.show()
#%%
fig, ax = plt.subplots(5, 6, figsize = (10, 8))
axes = ax.flatten()
for ind, key in enumerate(dict_seqname_to_stat.keys()):
    axes[ind].hist(dict_seqname_to_stat[key]["qv_by_read"], bins = np.linspace(0, 60, 51), color = "gray")
    if ind < 24:
        axes[ind].set_xticks([15])
    else:
        axes[ind].set_xticks([0, 15, 50])
        axes[ind].set_xlabel("QV per Read")
    if ind % 6 == 0:
        axes[ind].set_ylabel("# Reads")
    # axes[ind].set_title(key, fontsize = 10)
    axes[ind].ticklabel_format(axis = 'y', style = "sci", scilimits = (0, 0))
    axes[ind].axvline(15, linewidth = 1.5, color = "firebrick", linestyle = "--")
axes[-1].axis("off")
plt.tight_layout()
plt.show()
# %%
dict_seqname_to_ulqv = dict()
for key in dict_seqname_to_stat.keys():
    list_zip_length_to_qv = list(zip(dict_seqname_to_stat[key]["length_by_read"], dict_seqname_to_stat[key]["qv_by_read"]))
    list_zip_length_to_qv_only_ul = list(filter(lambda val: val[0] >= 100000, list_zip_length_to_qv))
    
    dict_seqname_to_ulqv[key] = list(map(lambda val: val[1], list_zip_length_to_qv_only_ul))

#%%
dict_seqname_to_ul_qvfiltered = dict()
for key in dict_seqname_to_stat.keys():
    list_zip_length_to_qv = list(zip(dict_seqname_to_stat[key]["length_by_read"], dict_seqname_to_stat[key]["qv_by_read"]))
    list_zip_length_to_qv_only_ul = list(filter(lambda val: val[0] >= 100000, list_zip_length_to_qv))
    
    list_zip_length_to_qv_only_ul_qvfilt = list(filter(lambda val: val[1] >= 15, list_zip_length_to_qv_only_ul))
    
    dict_seqname_to_ul_qvfiltered[key] = list_zip_length_to_qv_only_ul_qvfilt
#%%
fig, ax = plt.subplots(5, 6, figsize = (10, 8))
axes = ax.flatten()
for ind, key in enumerate(dict_seqname_to_ulqv.keys()):
    
    axes[ind].hist(dict_seqname_to_ulqv[key], bins = np.linspace(0, 60, 51), color = "gray")
    if ind < 24:
        axes[ind].set_xticks([15])
    else:
        axes[ind].set_xticks([0, 15, 50])
        axes[ind].set_xlabel("QV per Read")
    if ind % 6 == 0:
        axes[ind].set_ylabel("# Reads")
    # axes[ind].set_title(key, fontsize = 10)
    axes[ind].ticklabel_format(axis = 'y', style = "sci", scilimits = (0, 0))
    axes[ind].axvline(15, linewidth = 1.5, color = "firebrick", linestyle = "--")
axes[-1].axis("off")
plt.tight_layout()
plt.show()
#%%


plt.figure(figsize = (10, 3))
plt.bar(dict_seqname_to_stat.keys(),
        list(map(lambda key: sum(dict_seqname_to_stat[key]["length_by_read"]), dict_seqname_to_stat.keys())),
        color = "gray",
        label = "Total Bases",
        zorder = 0)

list_ul_bases = list(map(lambda key: 
                            sum(list(map(lambda val: val[0], dict_seqname_to_ul_qvfiltered[key])))
                            , dict_seqname_to_stat.keys()))

plt.bar(
    dict_seqname_to_stat.keys(),
    list_ul_bases,
    color = "firebrick",
    label = "UL Bases (Length ≥ 100kb & QV ≥ 15)",
    zorder = 3
)
for ind, key in enumerate(dict_seqname_to_stat.keys()):
    bases = sum(dict_seqname_to_stat[key]["length_by_read"])
    n_ul_bases = list_ul_bases[ind]
    plt.annotate(f"{n_ul_bases/bases*100:.2f} %", (ind, bases+1e9), rotation = 90, va = "bottom", ha = "center")
plt.xticks(rotation = -45, ha = "left")
plt.ylabel("Bases (bp)")
# plt.yscale("log")
plt.legend(loc = "center left", bbox_to_anchor = (1.01, 0.5))
plt.ylim(top = 1.4e11)
plt.show()

#%%
######################################
# Before vs After Trimming
######################################
#%%
from matplotlib import pyplot as plt

plt.rcParams["font.size"] = 14
plt.figure(figsize = (10, 3))

xticks = list(dict_seqname_to_stat.keys())

plt.bar(list(map(lambda val: val - 0.15, range(len(xticks)))),
        list(map(lambda key: np.mean(dict_seqname_to_stat[key]["qv_by_read"]), xticks)),
        color = "gray",
        width = 0.3,
        label = "Untrimmed")
plt.bar(list(map(lambda val: val + 0.15, range(len(xticks)))),
        list(map(lambda key: np.mean(dict_seqname_to_stat_trimmed[key]["qv_by_read"]), xticks)),
        color = "firebrick",
        width = 0.3, label = "Trimmed")


plt.xticks(range(len(xticks)), xticks, rotation = -45, ha = "left")
plt.ylabel("Mean QV")
plt.legend(loc = "center left", bbox_to_anchor = (1.01, 0.5))
plt.show()
# %%
plt.rcParams["font.size"] = 14
plt.figure(figsize = (10, 3))

xticks = list(dict_seqname_to_stat.keys())

plt.bar(list(map(lambda val: val - 0.15, range(len(xticks)))),
        list(map(lambda key: len(dict_seqname_to_stat[key]["length_by_read"]), xticks)),
        color = "gray",
        width = 0.3,
        label = "Untrimmed")
plt.bar(list(map(lambda val: val + 0.15, range(len(xticks)))),
        list(map(lambda key: len(dict_seqname_to_stat_trimmed[key]["length_by_read"]), xticks)),
        color = "firebrick",
        width = 0.3, label = "Trimmed")


plt.xticks(range(len(xticks)), xticks, rotation = -45, ha = "left")
plt.ylabel("# Reads")
plt.legend(loc = "center left", bbox_to_anchor = (1.01, 0.5))
plt.show()
# %%
plt.rcParams["font.size"] = 14
plt.figure(figsize = (10, 3))

xticks = list(dict_seqname_to_stat.keys())

plt.bar(list(map(lambda val: val - 0.15, range(len(xticks)))),
        list(map(lambda key: sum(dict_seqname_to_stat[key]["length_by_read"]), xticks)),
        color = "gray",
        width = 0.3,
        label = "Untrimmed")
plt.bar(list(map(lambda val: val + 0.15, range(len(xticks)))),
        list(map(lambda key: sum(dict_seqname_to_stat_trimmed[key]["length_by_read"]), xticks)),
        color = "firebrick",
        width = 0.3, label = "Trimmed")


plt.xticks(range(len(xticks)), xticks, rotation = -45, ha = "left")
plt.ylabel("# Base")
plt.legend(loc = "center left", bbox_to_anchor = (1.01, 0.5))
plt.show()

#%%
dict_read_to_qv_untrimmed = dict()
for key in dict_seqname_to_stat.keys():
    for readid in dict_seqname_to_stat[key]["readid"]:
        assert dict_read_to_qv_untrimmed.get(readid) == None, f"{key} > {readid}"
    dict_read_to_qv_untrimmed.update(
        dict(zip(dict_seqname_to_stat[key]["readid"], dict_seqname_to_stat[key]["qv_by_read"]))
    )
# %%
dict_read_to_qv_trimmed = dict()
for key in dict_seqname_to_stat_trimmed.keys():
    for readid in dict_seqname_to_stat_trimmed[key]["readid"]:
        assert dict_read_to_qv_trimmed.get(readid) == None, f"{key} > {readid}"
    dict_read_to_qv_trimmed.update(
        dict(zip(dict_seqname_to_stat_trimmed[key]["readid"], dict_seqname_to_stat_trimmed[key]["qv_by_read"]))
    )
# %%
import numpy as np

list_qv_diff = list(map(lambda key: dict_read_to_qv_trimmed[key] - dict_read_to_qv_untrimmed[key], dict_read_to_qv_trimmed.keys()))

plt.figure(figsize = (5, 3))
plt.hist(list_qv_diff, bins = 100, color = "gray")
plt.axvline(np.mean(list_qv_diff), linewidth = 1.5, color = "firebrick", linestyle = "--")
plt.xlabel("QV Difference")
plt.ylabel("# Reads")
plt.yscale("log")
plt.show()
#%%
plt.figure(figsize = (5, 3))
plt.scatter(
        list(map(lambda key: dict_read_to_qv_untrimmed[key], dict_read_to_qv_trimmed.keys())),
        list(map(lambda key: dict_read_to_qv_trimmed[key], dict_read_to_qv_trimmed.keys())),
        s = 5, color = "gray")
plt.plot([0, max(max(dict_read_to_qv_untrimmed.values()), max(dict_read_to_qv_trimmed.values()))], [0, max(max(dict_read_to_qv_untrimmed.values()), max(dict_read_to_qv_trimmed.values()))], linewidth = 1.5, color = "firebrick", linestyle = "--")
plt.xlabel("QV before Trimming")
plt.ylabel("QV after Trimming")
plt.show()
# %%
dict_seqname_to_ulqv_untrimmed = dict()
for key in dict_seqname_to_stat.keys():
    list_zip_length_to_qv = list(zip(dict_seqname_to_stat[key]["length_by_read"], dict_seqname_to_stat[key]["qv_by_read"]))
    list_zip_length_to_qv_only_ul = list(filter(lambda val: val[0] >= 100000, list_zip_length_to_qv))
    
    dict_seqname_to_ulqv_untrimmed[key] = list_zip_length_to_qv_only_ul

dict_seqname_to_ulqv_trimmed = dict()
for key in dict_seqname_to_stat_trimmed.keys():
    list_zip_length_to_qv = list(zip(dict_seqname_to_stat_trimmed[key]["length_by_read"], dict_seqname_to_stat_trimmed[key]["qv_by_read"]))
    list_zip_length_to_qv_only_ul = list(filter(lambda val: val[0] >= 100000, list_zip_length_to_qv))
    
    dict_seqname_to_ulqv_trimmed[key] = list_zip_length_to_qv_only_ul
#%%
plt.rcParams["font.size"] = 14
plt.figure(figsize = (10, 3))

xticks = list(dict_seqname_to_ulqv_untrimmed.keys())

plt.bar(list(map(lambda val: val - 0.15, range(len(xticks)))),
        list(map(lambda key: len(dict_seqname_to_ulqv_untrimmed[key]), xticks)),
        color = "gray",
        width = 0.3,
        label = "Untrimmed")
plt.bar(list(map(lambda val: val + 0.15, range(len(xticks)))),
        list(map(lambda key: len(dict_seqname_to_ulqv_trimmed[key]), xticks)),
        color = "firebrick",
        width = 0.3, label = "Trimmed")


plt.xticks(range(len(xticks)), xticks, rotation = -45, ha = "left")
plt.ylabel("# Reads")
plt.legend(loc = "center left", bbox_to_anchor = (1.01, 0.5))
plt.show()
# %%
plt.rcParams["font.size"] = 14
plt.figure(figsize = (10, 3))

xticks = list(dict_seqname_to_ulqv_untrimmed.keys())

plt.bar(list(map(lambda val: val - 0.15, range(len(xticks)))),
        list(map(lambda key: sum(list(map(lambda val: val[0], dict_seqname_to_ulqv_untrimmed[key]))), xticks)),
        color = "gray",
        width = 0.3,
        label = "Untrimmed")
plt.bar(list(map(lambda val: val + 0.15, range(len(xticks)))),
        list(map(lambda key: sum(list(map(lambda val: val[0], dict_seqname_to_ulqv_trimmed[key]))), xticks)),
        color = "firebrick",
        width = 0.3, label = "Trimmed")


plt.xticks(range(len(xticks)), xticks, rotation = -45, ha = "left")
plt.ylabel("# Base")
plt.legend(loc = "center left", bbox_to_anchor = (1.01, 0.5))
plt.show()
# %%
plt.rcParams["font.size"] = 14
plt.figure(figsize = (10, 3))

xticks = list(dict_seqname_to_ulqv_untrimmed.keys())

plt.bar(list(map(lambda val: val - 0.15, range(len(xticks)))),
        list(map(lambda key: np.mean(list(map(lambda val: val[1], dict_seqname_to_ulqv_untrimmed[key]))), xticks)),
        color = "gray",
        width = 0.3,
        label = "Untrimmed")
plt.bar(list(map(lambda val: val + 0.15, range(len(xticks)))),
        list(map(lambda key: np.mean(list(map(lambda val: val[1], dict_seqname_to_ulqv_trimmed[key]))), xticks)),
        color = "firebrick",
        width = 0.3, label = "Trimmed")


plt.xticks(range(len(xticks)), xticks, rotation = -45, ha = "left")
plt.ylabel("Mean QV")
plt.legend(loc = "center left", bbox_to_anchor = (1.01, 0.5))
plt.show()
# %%
dict_read_to_qv_ul_untrimmed = dict()
for key in dict_seqname_to_stat.keys():
    for readid, length, qv in zip(dict_seqname_to_stat[key]["readid"], dict_seqname_to_stat[key]["length_by_read"], dict_seqname_to_stat[key]["qv_by_read"]):
        if length >= 100000:
            dict_read_to_qv_ul_untrimmed[readid] = qv
# %%
dict_read_to_qv_ul_trimmed = dict()
for key in dict_seqname_to_stat_trimmed.keys():
    for readid, length, qv in zip(dict_seqname_to_stat_trimmed[key]["readid"], dict_seqname_to_stat_trimmed[key]["length_by_read"], dict_seqname_to_stat_trimmed[key]["qv_by_read"]):
        if length >= 100000:
            dict_read_to_qv_ul_trimmed[readid] = qv
# %%
list_qv_diff = list(map(lambda key: dict_read_to_qv_ul_trimmed[key] - dict_read_to_qv_ul_untrimmed[key], dict_read_to_qv_ul_trimmed.keys()))

plt.figure(figsize = (5, 3))
plt.hist(list_qv_diff, bins = 100, color = "gray")
plt.axvline(np.mean(list_qv_diff), linewidth = 1.5, color = "firebrick", linestyle = "--")
plt.xlabel("QV Difference")
plt.ylabel("# Reads")
plt.yscale("log")
plt.show()
# %%
plt.figure(figsize = (5, 3))
plt.scatter(
        list(map(lambda key: dict_read_to_qv_ul_untrimmed[key], dict_read_to_qv_ul_trimmed.keys())),
        list(map(lambda key: dict_read_to_qv_ul_trimmed[key], dict_read_to_qv_ul_trimmed.keys())),
        s = 5, color = "gray")
plt.plot([0, max(max(dict_read_to_qv_ul_untrimmed.values()), max(dict_read_to_qv_ul_trimmed.values()))], [0, max(max(dict_read_to_qv_ul_untrimmed.values()), max(dict_read_to_qv_ul_trimmed.values()))], linewidth = 1.5, color = "firebrick", linestyle = "--")
plt.xlabel("QV before Trimming")
plt.ylabel("QV after Trimming")
plt.ylim(10, 20)
plt.xlim(10, 20)
plt.show()
# %%
dict_seqname_to_ulqv_filtered_untrimmed = dict()
for key in dict_seqname_to_ulqv_untrimmed.keys():
    
    dict_seqname_to_ulqv_filtered_untrimmed[key] = list(filter(lambda val: val[1] >= 15, dict_seqname_to_ulqv_untrimmed[key]))


dict_seqname_to_ulqv_filtered_trimmed = dict()
for key in dict_seqname_to_ulqv_trimmed.keys():
    
    dict_seqname_to_ulqv_filtered_trimmed[key] = list(filter(lambda val: val[1] >= 15, dict_seqname_to_ulqv_trimmed[key]))
# %%
plt.rcParams["font.size"] = 14
plt.figure(figsize = (10, 3))

xticks = list(dict_seqname_to_ulqv_filtered_untrimmed.keys())

plt.bar(list(map(lambda val: val - 0.15, range(len(xticks)))),
        list(map(lambda key: len(dict_seqname_to_ulqv_filtered_untrimmed[key]), xticks)),
        color = "gray",
        width = 0.3,
        label = "Untrimmed")
plt.bar(list(map(lambda val: val + 0.15, range(len(xticks)))),
        list(map(lambda key: len(dict_seqname_to_ulqv_filtered_trimmed[key]), xticks)),
        color = "firebrick",
        width = 0.3, label = "Trimmed")


plt.xticks(range(len(xticks)), xticks, rotation = -45, ha = "left")
plt.ylabel("# Reads")
plt.legend(loc = "center left", bbox_to_anchor = (1.01, 0.5))
plt.show()
# %%
plt.rcParams["font.size"] = 14
plt.figure(figsize = (10, 3))

xticks = list(dict_seqname_to_ulqv_filtered_untrimmed.keys())

plt.bar(list(map(lambda val: val - 0.15, range(len(xticks)))),
        list(map(lambda key: np.mean(list(map(lambda val: val[1], dict_seqname_to_ulqv_filtered_untrimmed[key]))), xticks)),
        color = "gray",
        width = 0.3,
        label = "Untrimmed")
plt.bar(list(map(lambda val: val + 0.15, range(len(xticks)))),
        list(map(lambda key: np.mean(list(map(lambda val: val[1], dict_seqname_to_ulqv_filtered_trimmed[key]))), xticks)),
        color = "firebrick",
        width = 0.3, label = "Trimmed")


plt.xticks(range(len(xticks)), xticks, rotation = -45, ha = "left")
plt.ylabel("Mean QV")
plt.legend(loc = "center left", bbox_to_anchor = (1.01, 0.5))
plt.show()
# %%
from matplotlib import pyplot as plt
import numpy as np

fig, ax = plt.subplots(5, 6, figsize = (10, 8))
# plt.subplots_adjust(wspace = )
axes = ax.flatten()
for ind, key in enumerate(dict_seqname_to_stat.keys()):
    
    axes[ind].hist(dict_seqname_to_stat[key]["qv_by_read"], bins = np.linspace(0, 60, 51), color = "#2ca02c", alpha = 0.4, label = "Raw", density = True)
    # axes[ind].hist(list(map(lambda val: val[1], dict_seqname_to_ulqv_untrimmed[key])), bins = np.linspace(0, 60, 51), color = "#2ca02c", alpha = 0.4, label = "Untrim", density = True)
    
    axes[ind].hist(list(map(lambda val: val[1], dict_seqname_to_ulqv_untrimmed[key])), bins = np.linspace(0, 60, 51), color = "#9467bd", alpha = 0.4, label = "UL", density = True)
    # axes[ind].hist(list(map(lambda val: val[1], dict_seqname_to_ulqv_trimmed[key])), bins = np.linspace(0, 60, 51), color = "#9467bd", alpha = 0.4, label = "Trimmed", density = True)
    if ind < 24:
        axes[ind].set_xticks([])
    else:
        axes[ind].set_xticks([0, 15, 50])
        axes[ind].set_xlabel("QV")
    if ind % 6 == 0:
        axes[ind].set_ylabel("Ratio")
    axes[ind].set_yticks([])
    # axes[ind].set_title(key, fontsize = 10)
    axes[ind].ticklabel_format(axis = 'y', style = "sci", scilimits = (0, 0))
    axes[ind].axvline(15, linewidth = 1.5, color = "firebrick", linestyle = "--")
    
    if ind == 0:
        legend = axes[ind].legend()
        axes[-1].legend(handles=legend.legendHandles, labels=[t.get_text() for t in legend.texts], loc='center', bbox_to_anchor = (0.5, 0.5))
        axes[ind].get_legend().remove()

axes[-1].axis("off")
# plt.tight_layout()
plt.show()
# %%
qv_threshold = 15

dict_reads_qvthreshold_raw = dict()
for key in dict_seqname_to_stat.keys():
    list_qvs = dict_seqname_to_stat[key]["qv_by_read"]
    
    ratio_qvs_under_thresh = sum(list(map(lambda val: val < qv_threshold, list_qvs))) / len(list_qvs)
    ratio_qvs_over_thresh = sum(list(map(lambda val: val >= qv_threshold, list_qvs))) / len(list_qvs)
    
    dict_reads_qvthreshold_raw[key] = {
        "Under":ratio_qvs_under_thresh,
        "Over":ratio_qvs_over_thresh
    }

dict_reads_qvthreshold_ul = dict()
for key in dict_seqname_to_ulqv_untrimmed.keys():
    list_qvs = list(map(lambda val: val[1], dict_seqname_to_ulqv_untrimmed[key]))
    
    ratio_qvs_under_thresh = sum(list(map(lambda val: val < qv_threshold, list_qvs))) / len(list_qvs)
    ratio_qvs_over_thresh = sum(list(map(lambda val: val >= qv_threshold, list_qvs))) / len(list_qvs)
    
    dict_reads_qvthreshold_ul[key] = {
        "Under":ratio_qvs_under_thresh,
        "Over":ratio_qvs_over_thresh
    }
# %%
list_qv_under_raw = list(map(lambda key: dict_reads_qvthreshold_raw[key]["Under"], dict_reads_qvthreshold_raw.keys()))
list_qv_under_ul = list(map(lambda key: dict_reads_qvthreshold_ul[key]["Under"], dict_reads_qvthreshold_raw.keys()))

list_qv_over_raw = list(map(lambda key: dict_reads_qvthreshold_raw[key]["Over"], dict_reads_qvthreshold_raw.keys()))
list_qv_over_ul = list(map(lambda key: dict_reads_qvthreshold_ul[key]["Over"], dict_reads_qvthreshold_raw.keys()))
# %%
plt.figure(figsize = (1.5, 3))
plt.boxplot([list_qv_under_raw, list_qv_under_ul], zorder = 1, widths = 0.35, showfliers = False)
for ratio_raw, ratio_ul in zip(list_qv_under_raw, list_qv_under_ul):
    plt.plot([1, 2], [ratio_raw, ratio_ul], color = "gray")
plt.scatter([1]*len(list_qv_under_raw), list_qv_under_raw, zorder = 3, label = "Raw", color = "#2ca02c", s = 20)
plt.scatter([2]*len(list_qv_under_ul), list_qv_under_ul, zorder = 3, label = "UL", color = "#9467bd", s = 20)

plt.xticks([1, 2], ["Raw", "UL"])
plt.ylabel("Ratio")
plt.show()
# %%
plt.figure(figsize = (1.5, 3))
plt.boxplot([list_qv_over_raw, list_qv_over_ul], zorder = 1, widths = 0.35, showfliers = False)
for ratio_raw, ratio_ul in zip(list_qv_over_raw, list_qv_over_ul):
    plt.plot([1, 2], [ratio_raw, ratio_ul], color = "gray")
plt.scatter([1]*len(list_qv_over_raw), list_qv_over_raw, zorder = 3, label = "Raw", color = "#2ca02c", s = 20)
plt.scatter([2]*len(list_qv_over_ul), list_qv_over_ul, zorder = 3, label = "UL", color = "#9467bd", s = 20)

plt.xticks([1, 2], ["Raw", "UL"])
plt.ylabel("Ratio")
plt.show()
# %%
from matplotlib import pyplot as plt
import numpy as np

plt.rcParams["font.size"] = 14
plt.figure(figsize = (10, 3))

xticks = list(dict_seqname_to_stat.keys())

plt.bar(list(map(lambda val: val - 0.15, range(len(xticks)))),
        list(map(lambda key: sum(dict_seqname_to_stat[key]["n_LowQ_base"])/ sum(dict_seqname_to_stat[key]["length_by_read"])*100, xticks)),
        color = "gray",
        width = 0.3,
        label = "Untrimmed")
plt.bar(list(map(lambda val: val + 0.15, range(len(xticks)))),
        list(map(lambda key: sum(dict_seqname_to_stat_trimmed[key]["n_LowQ_base"])/ sum(dict_seqname_to_stat_trimmed[key]["length_by_read"])*100, xticks)),
        color = "firebrick",
        width = 0.3, label = "Trimmed")


plt.xticks(range(len(xticks)), xticks, rotation = -45, ha = "left")
plt.ylabel("<Q10 Base (%)")
plt.legend(loc = "center left", bbox_to_anchor = (1.01, 0.5))
plt.show()
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
        
    cnt_front_unt, cnt_end_unt = get_lowquality_base_count_from_end(dict_seqname_to_stat[seqname_repr], 5000)
    cnt_front_t, cnt_end_t = get_lowquality_base_count_from_end(dict_seqname_to_stat_trimmed[seqname_repr], 5000)
    
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
# %%
######################################
with open(f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/ONT/Adapter_Filtered/Stat_pickle/Edge_Lowq_Stat.Adapt_Filt.pickle", "rb") as fwb:
    dict_seqname_to_lowqcnt_untrimmed = pickle.load(fwb)
    
with open(f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/ONT/Edge_Trimmed/Stat_pickle/Edge_Lowq_Stat.Adapt_Filt.Edge_Trim_Front_1500_Tail_1500.pickle", "rb") as fwb:
    dict_seqname_to_lowqcnt_trimmed = pickle.load(fwb)
# %%
def get_length_count(list_length, check_length):
    count = np.zeros(check_length)
    for length in list_length:
        count += np.arange(check_length) < length
    return count

fig, ax = plt.subplots(1, 2, figsize = (12, 4))
for ind, (key, dict_lowq_cnt) in enumerate(dict_seqname_to_lowqcnt_trimmed.items()):
    cnt_front = dict_lowq_cnt["cnt_front"]
    cnt_end = dict_lowq_cnt["cnt_end"]
    check_length = 5000
    n_read = dict_lowq_cnt["n_read"]
    mean_perc_lowquality = dict_lowq_cnt["mean_perc_lowquality"]
    
    cnt_read_length = get_length_count(dict_seqname_to_stat_trimmed[key]["length_by_read"], check_length)
    
    ax[0].plot(range(check_length), cnt_front / cnt_read_length * 100, color = sns.color_palette("hls", len(dict_seqname_to_lowqcnt_untrimmed))[ind])
    # ax[0].axhline(mean_perc_lowquality, linewidth = 1.5, linestyle = '--', color = "firebrick")

    ax[1].plot(range(check_length), cnt_end / cnt_read_length * 100, color = sns.color_palette("hls", len(dict_seqname_to_lowqcnt_untrimmed))[ind])
    # ax[1].axhline(mean_perc_lowquality, linewidth = 1.5, linestyle = '--', color = "firebrick")
ax[0].set_xlabel("Coordinate (bp)")
ax[0].set_ylabel("Percentage Low Quality Base (%)")
ax[0].set_title("Front")
ax[1].set_xlabel("Coordinate (bp)")
ax[1].set_ylabel("Percentage Low Quality Base (%)")
ax[1].set_title("End")
plt.show()
# %%
