#%%
import sys,os,re
from copy import deepcopy

import pandas as pd
from pathlib import Path

sys.path.append(str(Path(__file__).parents[1]))
from fasta_parser import read_fasta_as_dict, save_dict_fasta_into_file


def run(path_ntlink_agp, path_telomere_end_info, path_fasta_from_ntlink, path_fasta_ntlink_input, path_fasta_teloonly, path_save_adjusted_fasta):
    os.makedirs(os.path.dirname(path_save_adjusted_fasta), exist_ok=True)

    table_agp = pd.read_csv(path_ntlink_agp, sep = '\t', names = ["obj", "obj_beg", "obj_end", "part_number", "comp_type", "comp_id", "comp_beg", "comp_end", "ori"])
    assert len(table_agp[~table_agp["comp_type"].isin(['U', 'N'])]["comp_id"].unique()) == table_agp[~table_agp["comp_type"].isin(['U', 'N'])].shape[0]

    table_telo_info = pd.read_csv(path_telomere_end_info, sep = '\t')

    table_agp_joined = table_agp[table_agp["obj"].apply(lambda val: val.startswith("ntLink"))]

    dict_ntlink_to_telo_ori = check_ntlink_scaffold_telomere_orientation(table_agp_joined, table_telo_info)
    dict_ntlink_to_teloclassify, dict_ntlink_to_unmatched = detect_ntlink_scaffolds_considering_telomere(dict_ntlink_to_telo_ori)
    
    dict_ntlink_to_teloclassify_corr, dict_ntlink_to_unmatched_corr = adjust_teloclassified_context_with_endgap(deepcopy(dict_ntlink_to_teloclassify), deepcopy(dict_ntlink_to_unmatched), table_agp)
    
    dict_wrong_t2t_link = check_teloboth_polarization(dict_ntlink_to_teloclassify_corr, table_agp, table_telo_info)
    print(f"Found {sum(list(map(len, dict_wrong_t2t_link.values())))} Unexpected T2T contexts (Ignore this linkage)")

    dict_ntlink = read_fasta_as_dict(path_fasta_from_ntlink)
    dict_scf_input = read_fasta_as_dict(path_fasta_ntlink_input)
    dict_scf_teloonly = read_fasta_as_dict(path_fasta_teloonly)

    dict_adjusted_scaffold = adjust_ntlink_scaffold_sequence_by_considering_telomere(dict_ntlink_to_teloclassify_corr, dict_ntlink_to_unmatched_corr, dict_wrong_t2t_link, table_agp, table_telo_info, dict_ntlink, dict_scf_input, dict_scf_teloonly)

    save_dict_fasta_into_file(dict_adjusted_scaffold, path_save_adjusted_fasta)

def check_ntlink_scaffold_telomere_orientation(table_agp, table_telo_info):
    dict_scaffold_to_teloend = dict(zip(table_telo_info["Scaffold"], table_telo_info["Telomere_End"]))

    dict_ntlink_to_telo_orientation = dict()

    for ntlink_id in table_agp["obj"].unique():
        table_agp_single_link = table_agp[table_agp["obj"] == ntlink_id].sort_values(by = ["obj_beg"])

        list_telo_ori = list()
        for _, row in table_agp_single_link.iterrows():
            comp_part_number = row["part_number"]
            telo_left = 0
            telo_right = 0
            if row["comp_type"] in ['U', 'N']:
                pass
            else:
                component_id = row["comp_id"]
                if "telomasked" in component_id:
                    joined_orientation = row["ori"]
                    component_id_wo_maskingname = component_id.replace("_telomasked", '')
                    teloend_ori = dict_scaffold_to_teloend[component_id_wo_maskingname]
                    if teloend_ori == "front":
                        if joined_orientation == '+':
                            telo_left = 1
                        else:
                            telo_right = 1
                    else:
                        if joined_orientation == '+':
                            telo_right = 1
                        else:
                            telo_left = 1
                else:
                    pass
            list_telo_ori.append(f"<{telo_left}-{comp_part_number}-{telo_right}>")
        dict_ntlink_to_telo_orientation[ntlink_id] = ''.join(list_telo_ori)
    return dict_ntlink_to_telo_orientation

def detect_ntlink_scaffolds_considering_telomere(dict_ntlink_to_telo_ori):
    re_telo_both = re.compile("<1-[0-9]+-0>(<0-[0-9]+-0>)*<0-[0-9]+-1>")
    re_telo_left = re.compile("<1-[0-9]+-0>(<0-[0-9]+-0>)+")
    re_telo_right = re.compile("(<0-[0-9]+-0>)+<0-[0-9]+-1>")
    re_notelo = re.compile("(<0-[0-9]+-0>)+")

    dict_ntlink_to_teloclassify = dict()
    dict_ntlink_to_unmatched = dict()
    
    for ntlink_id, telo_ori_info in dict_ntlink_to_telo_ori.items():
        teloinfo_tmp = deepcopy(telo_ori_info)

        dict_ntlink_to_teloclassify[ntlink_id] = {
            "Both":list(),
            "Left":list(),
            "Right":list(),
            "None":list()
        }
        for telo_both_found in re_telo_both.finditer(teloinfo_tmp):
            matched_context = telo_both_found.group(0)
            dict_ntlink_to_teloclassify[ntlink_id]["Both"].append(matched_context)
        for context in dict_ntlink_to_teloclassify[ntlink_id]["Both"]:
            teloinfo_tmp = teloinfo_tmp.replace(context, 'X')
            
        for telo_left_found in re_telo_left.finditer(teloinfo_tmp):
            matched_context = telo_left_found.group(0)
            dict_ntlink_to_teloclassify[ntlink_id]["Left"].append(matched_context)
        for context in dict_ntlink_to_teloclassify[ntlink_id]["Left"]:
            teloinfo_tmp = teloinfo_tmp.replace(context, 'X')
            
        for telo_right_found in re_telo_right.finditer(teloinfo_tmp):
            matched_context = telo_right_found.group(0)
            dict_ntlink_to_teloclassify[ntlink_id]["Right"].append(matched_context)
        for context in dict_ntlink_to_teloclassify[ntlink_id]["Right"]:
            teloinfo_tmp = teloinfo_tmp.replace(context, 'X')
            
        for telo_none_found in re_notelo.finditer(teloinfo_tmp):
            matched_context = telo_none_found.group(0)
            dict_ntlink_to_teloclassify[ntlink_id]["None"].append(matched_context)
        for context in dict_ntlink_to_teloclassify[ntlink_id]["None"]:
            teloinfo_tmp = teloinfo_tmp.replace(context, 'X')
        dict_ntlink_to_unmatched[ntlink_id] = teloinfo_tmp
    return dict_ntlink_to_teloclassify, dict_ntlink_to_unmatched

def adjust_teloclassified_context_with_endgap(dict_ntlink_to_teloclassify, dict_ntlink_to_unmatched, table_agp):
    re_anypart = re.compile("<[0,1]-[0-9]+-[0,1]>")
    
    for ntlink_id in dict_ntlink_to_teloclassify.keys():
        table_agp_thislink = table_agp[table_agp["obj"] == ntlink_id]
        list_parts_of_gap = table_agp_thislink[table_agp_thislink["comp_type"] == 'N']["part_number"].to_list()
        
        for classify_type, list_context in list(dict_ntlink_to_teloclassify[ntlink_id].items()):
            list_adjusted_context = list()
            for context in list_context:
                list_linked_component_number = get_linked_context_part_number_from_context_joined(context)
                list_context_part_per_component = re_anypart.findall(context)
                dict_component_number_to_context = dict(zip(list_linked_component_number, list_context_part_per_component))
                isgapcontext = ''.join(
                    list(map(
                        str, 
                        list(map(lambda compnumber: compnumber if compnumber not in list_parts_of_gap else 'N', list_linked_component_number))
                    ))
                )
                isgapcontext_without_endgap = isgapcontext.strip('N')
                if isgapcontext == isgapcontext_without_endgap:
                    list_adjusted_context.append(context)
                else:
                    print(f"Adjust End-gap Context: {ntlink_id}, {isgapcontext} -> {isgapcontext_without_endgap}")
                    if len(isgapcontext_without_endgap) == 1:
                        dict_ntlink_to_unmatched[ntlink_id] += dict_component_number_to_context[int(isgapcontext_without_endgap)]
                    elif len(isgapcontext_without_endgap) == 0:
                        pass
                    else:
                        list_adjusted_context.append(
                            ''.join(list(map(dict_component_number_to_context.__getitem__, list(map(int, isgapcontext_without_endgap)))))
                        )
            dict_ntlink_to_teloclassify[ntlink_id][classify_type] = list_adjusted_context
    return dict_ntlink_to_teloclassify, dict_ntlink_to_unmatched

def check_teloboth_polarization(dict_ntlink_to_teloclassify, table_agp, table_telo_info):
    dict_ori_to_polarize = {'-' : -1, '+' : 1}
    dict_sequencedirection_to_polarize = {"forward" : 1, "reverse" : -1}

    dict_seqname_to_telomeredirection = dict(zip(table_telo_info["Scaffold"], table_telo_info["Telomere_Sequence_Direction"]))

    dict_ntlink_wrong = dict()

    for ntlink_id in dict_ntlink_to_teloclassify.keys():
        table_agp_ntlinkid = table_agp[table_agp["obj"] == ntlink_id]
        for teloboth_context in dict_ntlink_to_teloclassify[ntlink_id]["Both"]:
            list_linked_component_number = get_linked_context_part_number_from_context_joined(teloboth_context)

            first_component = table_agp_ntlinkid[table_agp_ntlinkid["part_number"] == list_linked_component_number[0]]
            last_component = table_agp_ntlinkid[table_agp_ntlinkid["part_number"] == list_linked_component_number[-1]]

            assert first_component.shape[0] == 1, first_component
            assert last_component.shape[0] == 1, last_component

            dict_first_component = first_component.iloc[0, :].to_dict()
            dict_last_component = last_component.iloc[0, :].to_dict()

            first_component_original_name = dict_first_component["comp_id"].replace("_telomasked", '')
            first_component_telo_direction = dict_sequencedirection_to_polarize[dict_seqname_to_telomeredirection[first_component_original_name]] * dict_ori_to_polarize[dict_first_component["ori"]]

            last_component_original_name = dict_last_component["comp_id"].replace("_telomasked", '')
            last_component_telo_direction = dict_sequencedirection_to_polarize[dict_seqname_to_telomeredirection[last_component_original_name]] * dict_ori_to_polarize[dict_last_component["ori"]]

            if first_component_telo_direction + last_component_telo_direction != 0:
                if dict_ntlink_wrong.get(ntlink_id) == None:
                    dict_ntlink_wrong[ntlink_id] = list()
                dict_ntlink_wrong[ntlink_id].append(teloboth_context)
    return dict_ntlink_wrong

def adjust_ntlink_scaffold_sequence_by_considering_telomere(dict_ntlink_to_teloclassify, dict_ntlink_to_unmatched, dict_wrong_t2t_link, table_agp, table_telo_info, dict_ntlink, dict_scf_input, dict_scf_teloonly):
    dict_ntlink_adjusted = dict()

    set_ntlink_scfnames_in_agp = set(table_agp["obj"].unique())
    for scfname, scf_sequence in dict_ntlink.items():
        scfname = scfname.split()[0]
        assert scfname in set_ntlink_scfnames_in_agp, f"{scfname} not exists in .agp file"
        if dict_ntlink_to_teloclassify.get(scfname) == None:
            assert "ntLink" not in scfname, scfname

            original_scfname, original_scf_sequence = get_original_sequence_of_unntlinked_scaffold(table_agp, scfname, dict_scf_input, dict_scf_teloonly, table_telo_info)
            dict_ntlink_adjusted[original_scfname] = original_scf_sequence
        else:
            table_agp_of_thislink = table_agp[table_agp["obj"] == scfname]
            for ind_scf, both_telomere_context in enumerate(dict_ntlink_to_teloclassify[scfname]["Both"]):
                if both_telomere_context in dict_wrong_t2t_link.get(scfname, list()):
                    continue
                scfname_adjusted = f"{scfname}.Both.{ind_scf}"
                list_part_number = get_linked_context_part_number_from_context_joined(both_telomere_context)
                adjusted_scf_sequence = recover_trimmed_sequence_of_ntlink_with_part_numbers(table_agp_of_thislink, list_part_number, scf_sequence, dict_scf_input, n_base_check = 100)

                strand_first_component = table_agp_of_thislink[table_agp_of_thislink["part_number"] == list_part_number[0]].iloc[0, :].to_dict()["ori"]
                strand_last_component = table_agp_of_thislink[table_agp_of_thislink["part_number"] == list_part_number[-1]].iloc[0, :].to_dict()["ori"]
                telo_left = get_telomere_sequence_from_agp(table_agp_of_thislink, list_part_number[0], table_telo_info, dict_scf_teloonly, orient = strand_first_component)
                telo_right = get_telomere_sequence_from_agp(table_agp_of_thislink, list_part_number[-1], table_telo_info, dict_scf_teloonly, orient = strand_last_component)

                adjusted_scf_sequence = insert_telomere_sequence_in_range(adjusted_scf_sequence, telo_left, 0, len(telo_left))
                adjusted_scf_sequence = insert_telomere_sequence_in_range(adjusted_scf_sequence, telo_right, len(adjusted_scf_sequence)-len(telo_right), len(adjusted_scf_sequence))
                assert dict_ntlink_adjusted.get(scfname_adjusted) == None
                dict_ntlink_adjusted[scfname_adjusted] = adjusted_scf_sequence

            for ind_scf, left_telomere_context in enumerate(dict_ntlink_to_teloclassify[scfname]["Left"]):
                scfname_adjusted = f"{scfname}.Left.{ind_scf}"
                list_part_number = get_linked_context_part_number_from_context_joined(left_telomere_context)
                adjusted_scf_sequence = recover_trimmed_sequence_of_ntlink_with_part_numbers(table_agp_of_thislink, list_part_number, scf_sequence, dict_scf_input, n_base_check = 100)

                strand_first_component = table_agp_of_thislink[table_agp_of_thislink["part_number"] == list_part_number[0]].iloc[0, :].to_dict()["ori"]
                telo_left = get_telomere_sequence_from_agp(table_agp_of_thislink, list_part_number[0], table_telo_info, dict_scf_teloonly, orient = strand_first_component)

                adjusted_scf_sequence = insert_telomere_sequence_in_range(adjusted_scf_sequence, telo_left, 0, len(telo_left))
                assert dict_ntlink_adjusted.get(scfname_adjusted) == None
                dict_ntlink_adjusted[scfname_adjusted] = adjusted_scf_sequence
                
            for ind_scf, right_telomere_context in enumerate(dict_ntlink_to_teloclassify[scfname]["Right"]):
                scfname_adjusted = f"{scfname}.Right.{ind_scf}"
                list_part_number = get_linked_context_part_number_from_context_joined(right_telomere_context)
                adjusted_scf_sequence = recover_trimmed_sequence_of_ntlink_with_part_numbers(table_agp_of_thislink, list_part_number, scf_sequence, dict_scf_input, n_base_check = 100)
                
                strand_last_component = table_agp_of_thislink[table_agp_of_thislink["part_number"] == list_part_number[-1]].iloc[0, :].to_dict()["ori"]
                telo_right = get_telomere_sequence_from_agp(table_agp_of_thislink, list_part_number[-1], table_telo_info, dict_scf_teloonly, orient = strand_last_component)

                adjusted_scf_sequence = insert_telomere_sequence_in_range(adjusted_scf_sequence, telo_right, len(adjusted_scf_sequence)-len(telo_right), len(adjusted_scf_sequence))
                assert dict_ntlink_adjusted.get(scfname_adjusted) == None
                dict_ntlink_adjusted[scfname_adjusted] = adjusted_scf_sequence
                
            for ind_scf, none_telomere_context in enumerate(dict_ntlink_to_teloclassify[scfname]["None"]):
                scfname_adjusted = f"{scfname}.None.{ind_scf}"
                list_part_number = get_linked_context_part_number_from_context_joined(none_telomere_context)
                adjusted_scf_sequence = recover_trimmed_sequence_of_ntlink_with_part_numbers(table_agp_of_thislink, list_part_number, scf_sequence, dict_scf_input, n_base_check = 100)

                assert dict_ntlink_adjusted.get(scfname_adjusted) == None
                dict_ntlink_adjusted[scfname_adjusted] = adjusted_scf_sequence
                
    
    dict_original_fasta_of_unmatched_scaffolds = dict()
    for ntlink_id, context_unmatched in dict_ntlink_to_unmatched.items():
        context_unmatched_scfonly = context_unmatched.replace('X', '')
        if context_unmatched_scfonly:
            list_linked_component_number = get_linked_context_part_number_from_context_joined(context_unmatched_scfonly)
            dict_original_fasta_of_this_link = get_original_sequence_of_ntlink_with_part_numbers(table_agp, ntlink_id, list_linked_component_number, dict_scf_input, dict_scf_teloonly, table_telo_info)
            assert len(set(dict_original_fasta_of_unmatched_scaffolds.keys()) & set(dict_original_fasta_of_this_link)) == 0
            dict_original_fasta_of_unmatched_scaffolds.update(dict_original_fasta_of_this_link)
    assert len(set(dict_original_fasta_of_unmatched_scaffolds.keys()) & set(dict_ntlink_adjusted)) == 0
    dict_ntlink_adjusted.update(dict_original_fasta_of_unmatched_scaffolds)
    
    dict_original_fasta_of_wrong_link = dict()
    for wrong_scf_ntlink_id, list_wrong_link_context in dict_wrong_t2t_link.values():
        for wrong_link_context in list_wrong_link_context:
            list_linked_component_number = get_linked_context_part_number_from_context_joined(wrong_link_context)
            dict_original_fasta_of_this_link = get_original_sequence_of_ntlink_with_part_numbers(table_agp, wrong_scf_ntlink_id, list_linked_component_number, dict_scf_input, dict_scf_teloonly, table_telo_info)
            assert len(set(dict_original_fasta_of_wrong_link.keys()) & set(dict_original_fasta_of_this_link)) == 0
            dict_original_fasta_of_wrong_link.update(dict_original_fasta_of_this_link)
    assert len(set(dict_original_fasta_of_wrong_link.keys()) & set(dict_ntlink_adjusted)) == 0
    dict_ntlink_adjusted.update(dict_original_fasta_of_wrong_link)

    return dict_ntlink_adjusted
    
def get_linked_context_part_number_from_context_joined(context):
    list_linked_contexts = context.strip('<>').split('><')
    list_linked_component_number = list(map(lambda val: int(val.split('-')[1]), list_linked_contexts))
    return list_linked_component_number

def recover_trimmed_sequence_of_ntlink_with_part_numbers(table_agp_of_thislink, list_part_number, ntlink_sequence_full, dict_scf_input, n_base_check = 100):
    first_component = table_agp_of_thislink[table_agp_of_thislink["part_number"] == list_part_number[0]]
    last_component = table_agp_of_thislink[table_agp_of_thislink["part_number"] == list_part_number[-1]]

    assert first_component.shape[0] == 1, first_component
    assert last_component.shape[0] == 1, last_component

    dict_first_component = first_component.iloc[0, :].to_dict()
    dict_last_component = last_component.iloc[0, :].to_dict()

    ntlink_sequence = ntlink_sequence_full[dict_first_component["obj_beg"]-1:dict_last_component["obj_end"]]

    # First
    first_comp_name = dict_first_component["comp_id"]
    first_comp_beg = int(dict_first_component["comp_beg"])
    first_comp_end = int(dict_first_component["comp_end"])
    first_comp_sequence = dict_scf_input[first_comp_name]
    front_sequence_for_check = ntlink_sequence[:n_base_check]
    if dict_first_component["ori"] == '+':
        if first_comp_beg != 1:
            front_trimmed_sequence = first_comp_sequence[:first_comp_beg-1]
            front_sequence_check_from_input_scf = first_comp_sequence[first_comp_beg-1:first_comp_beg-1+n_base_check]
        else:
            front_trimmed_sequence = ''
            front_sequence_check_from_input_scf = first_comp_sequence[:n_base_check]

    elif dict_first_component["ori"] == '-':
        if first_comp_end != len(first_comp_sequence):
            front_trimmed_sequence = get_reverse_complement_sequence(first_comp_sequence[first_comp_end:])
            front_sequence_check_from_input_scf = get_reverse_complement_sequence(first_comp_sequence[first_comp_end-n_base_check:first_comp_end])
        else:
            front_trimmed_sequence = ''
            front_sequence_check_from_input_scf = get_reverse_complement_sequence(first_comp_sequence[-n_base_check:])
    else:
        raise Exception(f"Unexpected orientation: {dict_first_component}")
    assert front_sequence_check_from_input_scf == front_sequence_for_check

    # Last
    last_comp_name = dict_last_component["comp_id"]
    last_comp_beg = int(dict_last_component["comp_beg"])
    last_comp_end = int(dict_last_component["comp_end"])
    last_comp_sequence = dict_scf_input[last_comp_name]
    end_sequence_for_check = ntlink_sequence[-n_base_check:]
    if dict_last_component["ori"] == '+':
        if last_comp_end != len(last_comp_sequence):
            end_trimmed_sequence = last_comp_sequence[last_comp_end:]
            end_sequence_check_from_input_scf = last_comp_sequence[last_comp_end-n_base_check:last_comp_end]
        else:
            end_trimmed_sequence = ''
            end_sequence_check_from_input_scf = last_comp_sequence[-n_base_check:]

    elif dict_last_component["ori"] == '-':
        if last_comp_beg != 1:
            end_trimmed_sequence = get_reverse_complement_sequence(last_comp_sequence[:last_comp_beg-1])
            end_sequence_check_from_input_scf = get_reverse_complement_sequence(last_comp_sequence[last_comp_beg-1:last_comp_beg-1+n_base_check])
        else:
            end_trimmed_sequence = ''
            end_sequence_check_from_input_scf = get_reverse_complement_sequence(last_comp_sequence[:n_base_check])
    else:
        raise Exception(f"Unexpected orientation: {dict_last_component}")
    assert end_sequence_check_from_input_scf == end_sequence_for_check

    recovered_ntlink_sequence = front_trimmed_sequence+ntlink_sequence+end_trimmed_sequence
    return recovered_ntlink_sequence

def get_original_sequence_of_unntlinked_scaffold(table_agp, scaffold_name, dict_scf_input, dict_scf_teloonly, table_telo_info):
    dict_seqname_to_telostart = dict(zip(table_telo_info["Scaffold"], table_telo_info["Telomere_Masking_Start"]))
    dict_seqname_to_teloend = dict(zip(table_telo_info["Scaffold"], table_telo_info["Telomere_Masking_End"]))

    table_agp_scaffold = table_agp[table_agp["comp_id"] == scaffold_name]
    assert table_agp_scaffold.shape[0] == 1, scaffold_name
    dict_scaffold = table_agp_scaffold.iloc[0, :].to_dict()
    assert dict_scaffold["obj"] == dict_scaffold["comp_id"], scaffold_name

    sequence = dict_scf_input[scaffold_name]
    if scaffold_name.endswith("_telomasked"):
        scfname_original = scaffold_name.replace("_telomasked", '')
        telo_sequence = get_telomere_sequence_from_agp(table_agp_scaffold, 1, table_telo_info, dict_scf_teloonly)
        telo_start = dict_seqname_to_telostart[scfname_original]
        telo_end = dict_seqname_to_teloend[scfname_original]
        sequence = insert_telomere_sequence_in_range(dict_scf_input[scaffold_name], telo_sequence, telo_start, telo_end)

        scaffold_name = scfname_original
    return scaffold_name, sequence

def get_original_sequence_of_ntlink_with_part_numbers(table_agp, ntlink_id, list_part_number, dict_scf_input, dict_scf_teloonly, table_telo_info):
    dict_original_sequence = dict()

    dict_seqname_to_telostart = dict(zip(table_telo_info["Scaffold"], table_telo_info["Telomere_Masking_Start"]))
    dict_seqname_to_teloend = dict(zip(table_telo_info["Scaffold"], table_telo_info["Telomere_Masking_End"]))

    table_agp_id = table_agp[table_agp["obj"] == ntlink_id]
    table_agp_id_part = table_agp_id[table_agp_id["part_number"].isin(list_part_number)]
    table_agp_id_part_scf = table_agp_id_part[~table_agp_id_part["comp_type"].isin(['U', 'N'])]

    for scfname, part_number in zip(table_agp_id_part_scf["comp_id"], table_agp_id_part_scf["part_number"]):
        sequence_input = dict_scf_input[scfname]
        if scfname.endswith("telomasked"):
            scfname_original = scfname.replace("_telomasked", '')
            telo_sequence = get_telomere_sequence_from_agp(table_agp_id_part_scf, part_number, table_telo_info, dict_scf_teloonly)
            telo_start = dict_seqname_to_telostart[scfname_original]
            telo_end = dict_seqname_to_teloend[scfname_original]
            dict_original_sequence[scfname_original] = insert_telomere_sequence_in_range(sequence_input, telo_sequence, telo_start, telo_end)
        else:
            dict_original_sequence[scfname] = sequence_input
    return dict_original_sequence

def get_telomere_sequence_from_agp(table_agp, part_number, table_telo_info, dict_scf_telonly, orient = '+'):
    dict_seqname_to_teloonly_seqname = dict(zip(table_telo_info["Scaffold"], table_telo_info.apply(lambda row: f"{row['Scaffold']}_teloonly_{row['Telomere_Masking_Start']}_{row['Telomere_Masking_End']}", axis = 1)))

    table_agp_of_specific_component = table_agp[table_agp["part_number"] == part_number]
    assert table_agp_of_specific_component.shape[0] == 1, table_agp_of_specific_component
    dict_component = table_agp_of_specific_component.iloc[0, :].to_dict()
    assert dict_component["comp_id"].endswith("_telomasked"), dict_component

    compname_original = dict_component["comp_id"].replace("_telomasked", '')
    teloname = dict_seqname_to_teloonly_seqname[compname_original]

    telo_sequence = dict_scf_telonly[teloname]
    if orient == '-':
        telo_sequence = get_reverse_complement_sequence(telo_sequence)
    return telo_sequence

def insert_telomere_sequence_in_range(sequence, telo_sequence, range_start, range_end):
    assert set(sequence[range_start:range_end]) == {'N'}, "Seems this region is not masked!"
    if len(sequence) == range_end:
        sequence_recovered = sequence[:range_start] + telo_sequence
    elif range_start == 0:
        sequence_recovered = telo_sequence + sequence[range_end:]
    return sequence_recovered

def get_reverse_complement_sequence(sequence):
    dict_reverse = {
        'A' : 'T',
        'T' : 'A',
        'G' : 'C',
        'C' : 'G',
        'N' : 'N'
    }

    reversed_sequence = sequence[::-1]
    reverse_complement_sequence = list(map(dict_reverse.__getitem__, reversed_sequence))
    return ''.join(reverse_complement_sequence)
    

# %%

if __name__ == "__main__":
    hap = 2

    dict_hap_to_side = {
        1 : "p",
        2 : "m"
    }

    # path_ntlink_agp = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.asm.dip.hap{hap}.no_telo.telo_single_masked.ntLink.ONT_triobin_{dict_hap_to_side[hap]}.L50kb_QV20_Q30Perc80.sensitive_gapfilling.trimmed_scafs.agp"
    # path_telomere_end_info = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/Telomere_Detection/Fasta_Split_by_Telomere/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.asm.dip.hap{hap}.p_ctg.telo_single.telo_info.tsv"
    # path_fasta_from_ntlink = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.asm.dip.hap{hap}.p_ctg.no_telo.telo_single_masked.fa.k32.w100.z1000.ntLink.scaffolds.fa"
    # path_fasta_ntlink_input = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/Telomere_Detection/Fasta_Split_by_Telomere/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.asm.dip.hap{hap}.p_ctg.no_telo.telo_single_masked.fasta"
    # path_fasta_teloonly = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/Telomere_Detection/Fasta_Split_by_Telomere/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.asm.dip.hap{hap}.p_ctg.telo_single.masked_telo_only.fasta"
    # path_save_adjusted_fasta = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling/Adjusted/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.asm.dip.hap{hap}.p_ctg.no_telo.telo_single_masked.fa.k32.w100.z1000.ntLink.scaffolds.adjusted_telomere.fa"

    path_ntlink_agp = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.no_telo.telo_single_masked.ntLink.ONT_triobin_{dict_hap_to_side[hap]}.L50kb_QV20_Q30Perc80.sensitive_gapfilling.trimmed_scafs.agp"
    path_telomere_end_info = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/Telomere_Detection/Fasta_Split_by_Telomere/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg.telo_single.telo_info.tsv"
    path_fasta_from_ntlink = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg.no_telo.telo_single_masked.fa.k32.w100.z1000.ntLink.scaffolds.fa"
    path_fasta_ntlink_input = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/Telomere_Detection/Fasta_Split_by_Telomere/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg.no_telo.telo_single_masked.fasta"
    path_fasta_teloonly = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/Telomere_Detection/Fasta_Split_by_Telomere/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg.telo_single.masked_telo_only.fasta"
    path_save_adjusted_fasta = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling/Adjusted/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg.no_telo.telo_single_masked.fa.k32.w100.z1000.ntLink.scaffolds.adjusted_telomere.fa"

    run(path_ntlink_agp, path_telomere_end_info, path_fasta_from_ntlink, path_fasta_ntlink_input, path_fasta_teloonly, path_save_adjusted_fasta)
# %%
