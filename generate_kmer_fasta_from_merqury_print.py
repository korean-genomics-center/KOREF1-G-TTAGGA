#%%
from zlib_ng import gzip_ng
import os


def convert_merqury_kmer_print_into_fasta(path_merqury_print, path_save_fasta):
    os.makedirs(os.path.dirname(path_save_fasta), exist_ok = True)

    file_reader = open(path_merqury_print, 'r')
    file_writer = gzip_ng.open(path_save_fasta, 'wb')

    for line in file_reader:
        sequence, count = line.strip('\n').split()
        count = int(count)

        for ind in range(count):
            seqname = f"kmer_{sequence}.iter_{ind}"
            file_writer.write(f">{seqname}\n".encode())
            file_writer.write(f"{sequence}\n".encode())
    
    file_reader.close()
    file_writer.close()       


#%%
if __name__ == "__main__":
    from joblib import Parallel, delayed

    with Parallel(n_jobs = 2) as parallel:
        parallel(delayed(convert_merqury_kmer_print_into_fasta)(
            path_merqury_print = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/Binning_Compare/{parent}_Trimmed_F15_T5.noNbase.merged.k31.hapmer.kmer_count.txt",
            path_save_fasta = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/Binning_Compare/Fasta_from_Meryl_Print/{parent}_Trimmed_F15_T5.noNbase.merged.k31.hapmer.kmer_count.fq.gz"
        ) for parent in ["MATERNAL", "PATERNAL"])

# %%
