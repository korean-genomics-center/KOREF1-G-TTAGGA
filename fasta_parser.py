#%%

def read_fasta_as_dict(path_fasta):
    dict_fasta = dict()

    with open(path_fasta, 'r') as fr:
        while 1:
            header = fr.readline()
            if not header:
                break
            assert header.startswith('>')
            seqname = header[1:].strip('\n')
            sequence = fr.readline().strip('\n')

            assert dict_fasta.get(seqname) == None
            dict_fasta[seqname] = sequence
    return dict_fasta

def save_specific_seqnames_from_fasta_dict(dict_fasta, list_names, path_save):
    import os
    os.makedirs(os.path.dirname(path_save), exist_ok = True)
    with open(path_save, 'w') as fw:
        for name in list_names:
            fw.write(f">{name}\n")
            fw.write(dict_fasta[name] + '\n')
# %%
