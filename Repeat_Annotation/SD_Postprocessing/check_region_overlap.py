#%%
import os


def run(path_bedpe, path_save):
    with open(path_bedpe, 'r') as fr, open(path_save, 'w') as fw:
        for line in fr:
            values = line.strip('\n').split('\t')
            chr1, start1, end1, chr2, start2, end2 = values[:6]

            start1 = int(start1)
            end1 = int(end1)
            start2 = int(start2)
            end2 = int(end2)

            if chr1 == chr2:
                int_start = max(start1, start2)
                int_end = min(end1, end2)
                overlap = max(0, int_end - int_start)
            else:
                overlap = 0
            
            values.append(overlap)
            fw.write('\t'.join(list(map(str, values)))+'\n')

if __name__ == "__main__":
    for hap in [1, 2]:

        run(
            f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Repeat_Annotation/Segmental_Duplication/Biser/Hap{hap}/FinalASM.KOREF_S1.Trio_Specific.hap{hap}.Reordered.masked.SD",
            f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Repeat_Annotation/Segmental_Duplication/Biser/Hap{hap}/FinalASM.KOREF_S1.Trio_Specific.hap{hap}.Reordered.masked.SD.check_overlap"
        )