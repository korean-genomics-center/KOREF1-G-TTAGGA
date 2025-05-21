#%%
from collections import Counter
import gzip, subprocess
import numpy as np
from zlib_ng import gzip_ng

#%%

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

def trim_fastq_front_and_tail(path_fasta, path_save, is_gzip, trim_front, trim_tail, min_length):
    import numpy as np
    
    filereader = FASTA_READER(path_fasta, is_gzip)
    filewriter = FASTA_WRITER(path_save, is_gzip)
    
    set_atgc = set("ATGCN")
    
    n_passed = 0
    n_filtered_by_length = 0
    
    while 1:
        namevalue = filereader.readline()
        if namevalue == '':
            break
        if namevalue.startswith('@'):
            name = namevalue.strip("@\n")
        else:
            continue
        sequencevalue = filereader.readline()
        tmp = filereader.readline()
        qualityvalue = filereader.readline()
        
        sequence = sequencevalue.strip('\n').upper()
        quality = qualityvalue.strip('\n')
        
        assert set(sequence) <= set_atgc, name
        
        len_read = len(sequence)
        
        write_line = True
        if len_read < (trim_front + trim_tail + min_length):
            write_line = False
            n_filtered_by_length += 1
        if write_line:
            n_passed += 1
            filewriter.write(namevalue)
            filewriter.write(sequence[trim_front:-trim_tail]+'\n')
            filewriter.write(tmp)
            filewriter.write(quality[trim_front:-trim_tail]+'\n')
        
        n_processed = n_passed+n_filtered_by_length
        if (n_processed) % 10000 == 0:
            print(n_processed)
        
    filereader.close()
    filewriter.close()
    return (n_passed, n_filtered_by_length)

def run(path_fasta, path_save, is_gzip, trim_front, trim_tail, min_length):
    import os
    import numpy as np    
    
    n_passed, n_filtered_by_length = trim_fastq_front_and_tail(
        path_fasta, path_save, is_gzip, trim_front, trim_tail, min_length
    )
    
    with open(f"{path_save}.stat", 'w') as fw:
        fw.write(f"#Total : {n_passed+n_filtered_by_length}\n")
        fw.write(f"#Passed : {n_passed}\n")
        fw.write(f"#Filtered_by_Length : {n_filtered_by_length}\n")

    print("Done:", path_fasta, sep = ' ')
    
#%%


if __name__ == "__main__":
    list_seq_name = [
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
    trim_front = 1500
    trim_tail = 1500
    min_length = 50
    
    from joblib import Parallel, delayed
    
    with Parallel(n_jobs = len(list_seq_name)) as parallel:
        parallel(
            delayed(run)(
                f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/ONT/Adapter_Filtered/{seqname}.Adapt_Filt.fastq",
                f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/ONT/Edge_Trimmed/{seqname}.Adapt_Filt.Edge_Trim_Front_{trim_front}_Tail_{trim_tail}.fastq",
                False,
                trim_front,
                trim_tail,
                min_length
            ) for seqname in list_seq_name
        )