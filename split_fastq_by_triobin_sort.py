#%%
from zlib_ng import gzip_ng
import os
import pandas as pd
#%%

def main(path_fastq, path_triobin, path_output_prefix):
    os.makedirs(os.path.dirname(path_output_prefix), exist_ok = True)

    table_triobin = pd.read_csv(path_triobin, sep = '\t', names = "ctg err strongMixed sPat sMat weakMixed wPat1 wMat1 wPat2 wMat2".split())

    assert len(table_triobin["ctg"].unique()) == table_triobin.shape[0]

    dict_seqname_to_triobin = dict(zip(table_triobin["ctg"], table_triobin["err"].astype(str)))

    path_save_m = f"{path_output_prefix}.triobin_m.fq.gz"
    path_save_p = f"{path_output_prefix}.triobin_p.fq.gz"
    path_save_a = f"{path_output_prefix}.triobin_a.fq.gz"
    path_save_0 = f"{path_output_prefix}.triobin_0.fq.gz"


    is_gzip = path_fastq.endswith(".gz")
    if is_gzip:
        file_reader = gzip_ng.open(path_fastq, "rb")
    else:
        file_reader = open(path_fastq, 'r')

    dict_filewriter = {
        'm' : gzip_ng.open(path_save_m, "wb"),
        'p' : gzip_ng.open(path_save_p, "wb"),
        'a' : gzip_ng.open(path_save_a, "wb"),
        '0' : gzip_ng.open(path_save_0, "wb")
    }
    
    while 1:
        seqname = file_reader.readline()
        if not seqname:
            break
        if is_gzip:
            seqname_val = seqname.decode().strip('\n')
        else:
            seqname_val = seqname.strip('\n')
        assert seqname_val.startswith('@')

        sequence = file_reader.readline()
        tmp = file_reader.readline()
        quality = file_reader.readline()

        seqname_id = seqname_val[1:].split()[0]

        triobin_class = dict_seqname_to_triobin[seqname_id]

        if is_gzip:
            dict_filewriter[triobin_class].write(seqname)
            dict_filewriter[triobin_class].write(sequence)
            dict_filewriter[triobin_class].write(tmp)
            dict_filewriter[triobin_class].write(quality)
        else:
            dict_filewriter[triobin_class].write(seqname.encode())
            dict_filewriter[triobin_class].write(sequence.encode())
            dict_filewriter[triobin_class].write(tmp.encode())
            dict_filewriter[triobin_class].write(quality.encode())
    
    file_reader.close()
    for writer in dict_filewriter.values():
        writer.close()


if __name__ == "__main__":
    main(
        path_fastq = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/ONT/Merged/ONT_UL.Adapt_Filt.Edge_Trim_Front_1500_Tail_1500.Filter_L10kb_QV15.merged.fastq", 
        path_triobin = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/ONT/Parental_Split/ONT_UL.Adapt_Filt.Edge_Trim_Front_1500_Tail_1500.Filter_L10kb_QV15.merged.split_info.txt", 
        path_output_prefix = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/ONT/Parental_Split/ONT_UL.Adapt_Filt.Edge_Trim_Front_1500_Tail_1500.Filter_L10kb_QV15.merged.triobin_split")