# CEG
# haochenli
# Python3
# v0.2.5
import sys
import os
import pandas as pd
import re
from collections import Counter
import argparse
import multiprocessing


def read_table(file_name):
    # Read the file into a list
    out_table = []
    in_put = open(file_name, "r")
    for line1 in in_put:
        line1_rs1 = line1.rstrip("\n")
        line1_rs2 = line1_rs1.strip(" ")
        out_table.append(line1_rs2)
    in_put.close()
    return out_table


def read_dic(file_name):
    # Read a dictionary from a file
    out_dic = {}
    keys = []
    in_put = open(file_name, "r")
    for line2 in in_put:
        line2_sp = line2.split("\t")
        key = line2_sp[0]
        value = "\t".join(line2_sp[1:])
        out_dic[key] = value
        keys.append(key)
    in_put.close()
    return out_dic


def complete_path(input_path):
    # Ensure that the input is an absolute path
    if not os.path.isabs(input_path):
        input_path = os.path.join(os.getcwd(), input_path)
    return input_path


def enumerate_count(file_name):
    # Progress statistics
    with open(file_name) as f:
        for count, _ in enumerate(f, 1):
            pass
    return count


def get_file_matched_num(rex_str, net_f_table, output_matched_table_path):
    # Calculate the number matched
    with open(output_matched_table_path, 'a+', encoding='utf-8') as f:
        rex_obj = re.compile(rex_str, re.VERBOSE)
        count = 0
        for i in net_f_table:
            if rex_obj.findall(i):
                out_str = rex_str + "\t" + i
                count = count + 1
                print(out_str.rstrip("\n"), file=f)
        return count


def fill_dataframe_num_sp(loc_info, dataframe):
    # Fill in the matrix based on the count results
    for l1 in loc_info:
        l1_sp = l1.split("\t")
        row_loc = l1_sp[0]
        col_loc = l1_sp[1]
        value = l1_sp[2]
        dataframe.loc[row_loc, col_loc] = int(value)


def get_loc_table(count_table, sp_cls_table):
    # Obtain the filled information
    all_loc_table = []
    cls_table = read_table(sp_cls_table)
    cls_dic = {}
    for l1 in cls_table:
        if not l1.startswith("#"):
            l1_sp = l1.split("\t")
            cls_dic[l1_sp[0]] = l1_sp[5]
    for l2 in count_table:
        if "_" in l2:
            l2_sp1 = l2.split("\t")
            count = l2_sp1[1]
            l2_sp2 = l2_sp1[0].split("_")
            str1 = cls_dic.get(l2_sp2[0])
            str2 = cls_dic.get(l2_sp2[1])
            out_str_table = [str1, str2, count]
            all_loc_table.append(out_str_table)
        else:
            l2_sp1 = l2.split("\t")
            count = l2_sp1[1]
            str1 = cls_dic.get(l2_sp1[0])
            out_str_table = [str1, str1, count]
            all_loc_table.append(out_str_table)
    return all_loc_table


def get_flank_gene_in_bed(input_gene_id, input_bed_path, input_flanking_gene_length, output_flanking_gene_list_path):
    # Extract flanking genes ID from the bed file
    bed_table = read_table(input_bed_path)
    col1_table = []
    output_df = []
    for l1 in bed_table:
        l1_sp = l1.split("\t")
        if l1_sp[1] == input_gene_id:
            col1_table.append(l1_sp[0])

    bed_sp_table = []
    for l2 in bed_table:
        l2_sp = l2.split("\t")
        if l2_sp[0] == col1_table[0]:
            line_table = []
            for l2_2 in l2_sp:
                line_table.append(l2_2)
            bed_sp_table.append(line_table)

    df = pd.DataFrame(bed_sp_table, columns=['col1', 'col2', 'col3', 'col4'])
    df['col3'] = df['col3'].astype(int)
    sorted_df = df.sort_values(by='col3', ascending=True)
    sorted_df.reset_index(drop=True, inplace=True)
    print(sorted_df)
    match_index = sorted_df[sorted_df['col2'] == input_gene_id].index
    print(match_index)
    if match_index.empty:
        print(f"no math ID found: {input_gene_id}")
    else:
        for index in match_index:
            start = max(index - int(input_flanking_gene_length), 0)
            end = min(index + int(input_flanking_gene_length) + 1, len(sorted_df))
            output_df = sorted_df.iloc[start:end]
    output_df['col2'].to_csv(output_flanking_gene_list_path, index=False, header=False)


def stat_flanking_gene_syn_match_num(output_flanking_gene_list_path, input_syn_file, tem_output_path,
                                     output_matched_table):
    # Count the number of hits for flanking genes in the overall synteny network
    flanking_gene_list_table = read_table(output_flanking_gene_list_path)
    task_volume = enumerate_count(output_flanking_gene_list_path)
    syn_f_table = read_table(input_syn_file)
    with open(tem_output_path, 'a+', encoding='utf-8') as f1:
        i = 0
        for l3 in flanking_gene_list_table:
            i = i + 1
            print('find' + ' ' + str(i) + ' of ' + str(task_volume), end="\n", flush=True)
            out_str = l3 + "\t" + str(get_file_matched_num(l3, syn_f_table, output_matched_table))
            print('out str =' + out_str)
            print(out_str.rstrip("\n"), file=f1)


def stat_flanking_gene_syn_clade_match_num(output_flanking_gene_list_path, output_matched_table_path, sp_dic_path,
                                           output_sp_count_stat_path, output_sp_count_matrix_path):
    # Count the statistical results according to the species' classification information
    flanking_gene_list_table = read_table(output_flanking_gene_list_path)
    matched_table = read_table(output_matched_table_path)
    sp_dic = read_dic(sp_dic_path)
    for l4 in flanking_gene_list_table:
        print(f"Count the collinearity of {l4}")
        sp_abb_table = []
        clade_table = []
        replace_table = []
        for l5 in matched_table:
            l5_sp = l5.split("\t")
            if l5_sp[0] in l4:
                if l5_sp[0][0:4][-1] == '_':
                    replace_table.append(l5_sp[0][0:3])
                if sp_dic.get(l5_sp[0][0:3]) is not None:
                    replace_table.append(l5_sp[0][0:3])
                if sp_dic.get(l5_sp[0][0:4]) is not None:
                    replace_table.append(l5_sp[0][0:4])
                sp_abb_prefix1 = re.sub('-.*', '', l5_sp[1])
                sp_abb_prefix2 = re.sub(r'\d+', '', sp_abb_prefix1)
                sp_abb = sp_abb_prefix2.replace(replace_table[-1], '')
                print(sp_abb)
                if len(sp_abb) == 5:
                    sp_abb = sp_abb.replace('_', '')
                    sp_abb_table.append(sp_abb)
                if len(sp_abb) == 4:
                    if '_' in sp_abb:
                        sp_abb = sp_abb.replace('_', '')
                        sp_abb_table.append(sp_abb)
                    else:
                        sp_abb_table.append(sp_abb)
                if replace_table[-1] == sp_abb:
                    sp_abb_table.append(sp_abb)
        print(f'species_abbreviation:{sp_abb_table}')
        for l6 in sp_abb_table:
            try:
                sp_value = sp_dic.get(l6).split("\t")
                clade = sp_value[1]
                clade_table.append(clade)
            except AttributeError:
                pass
        count = Counter(clade_table)
        print(count)

        with open(output_sp_count_stat_path, 'a+', encoding='utf-8') as f2:
            for item, count in count.items():
                print(f"{l4}\t{item}\t{count}", file=f2)

    data_frame_row = read_table(output_flanking_gene_list_path)
    loc_table = read_table(output_sp_count_stat_path)
    # The content in the columns can be customized according to requirements.
    df_num = pd.DataFrame(columns=["Super-Rosids", "Super-Asterids", "Basal-Eudicots", "Monocots", "Magnoliids",
                                   'Basal-Angiosperm'], index=data_frame_row)
    df_num = df_num.fillna(0).infer_objects(copy=False)
    fill_dataframe_num_sp(loc_table, df_num)

    with open(output_sp_count_matrix_path, 'a+', encoding='utf-8') as f3:
        # When changing the columns, be sure to also update the header accordingly.
        header_str = f"ID\tSuper_Rosids\tSuper_Asterids\tBasal_Eudicots\tMonocots\tMagnoliids\tBasal_Angiosperm"
        print(header_str.rstrip("\n"), file=f3)
    df_num.to_csv(output_sp_count_matrix_path, mode='a', sep='\t', index=True, header=False)


def synteny_mapping(gene_id, input_bed, output_path, flanking_gene_length, input_syn_dic, input_sp_dic):
    # Implementation of the synteny_mapping function
    try:
        output_flanking_gene_list = f'{output_path}/{gene_id}_{flanking_gene_length}_flanking_gene.namelist'
        tem_output = f"{output_path}/{gene_id}_{flanking_gene_length}.tem2"
        output_matched_table = f"{output_path}/{gene_id}_{flanking_gene_length}.matched.tsv"
        output_sp_count_stat = f"{output_path}/{gene_id}_{flanking_gene_length}.stat.tsv"
        output_sp_count_matrix = f"{output_path}/{gene_id}_{flanking_gene_length}.stat.matrix.tsv"
        get_flank_gene_in_bed(gene_id, input_bed, flanking_gene_length, output_flanking_gene_list)
        stat_flanking_gene_syn_match_num(output_flanking_gene_list, input_syn_dic, tem_output, output_matched_table)
        stat_flanking_gene_syn_clade_match_num(output_flanking_gene_list, output_matched_table, input_sp_dic,
                                               output_sp_count_stat, output_sp_count_matrix)
    except IndexError:
        print(f'{gene_id} not found')
        pass


def synteny_mapping_mutiprocess(gene_table, bed, out_dir, size, input_syn_dic, species_namelist_file, num_processes):
    # Assign processes to different fragment search tasks
    pool = multiprocessing.Pool(processes=num_processes)
    pool.starmap(synteny_mapping, [(l1, bed, out_dir, size, input_syn_dic, species_namelist_file) for l1 in gene_table])
    pool.close()
    pool.join()


def main():
    text1 = '''Multi-species synteny mapping of chromosomal regions
            '''
    text2 = '''-i: list of species names; 
        input a list of file names for each species' pep and bed files
        '''
    text3 = '''-l: list of gene ID; input a list of gene ID;
        '''
    text4 = '''--bed: input target species' bed file; 
        '''
    text5 = '''-n: input Synnet_file; 
    input Synnet_file which generated through SynBuild
    '''
    text6 = '''-o: output_dir; 
        Output files to the specified directory
        '''
    text7 = '''-S: determine the block size (default=100); Determine the maximum number of genes to be retained upstream 
        and downstream of the target gene in the chromosome/scaffold.
        '''
    text8 = '''-p: number of threads: Number of threads for synteny mapping (default=1);
        only takes effect when there are multiple IDs in the input list.
        '''
    parser = argparse.ArgumentParser(description=text1)
    parser.add_argument('-i', '--species_namelist_file', type=str, required=True, help=text2)
    parser.add_argument('-l', '--id_list', type=str, required=True, help=text3)
    parser.add_argument('--bed', type=str, required=True, help=text4)
    parser.add_argument('-n', '--SynNet_f', type=str, required=True, help=text5)
    parser.add_argument('-o', '--out_dir', type=str, required=True, help=text6)
    parser.add_argument('-S', '--size', type=int, required=True, default=50, help=text7)
    parser.add_argument('-p', '--threads', type=int, default=1, help=text8)
    args = parser.parse_args()
    para = " ".join(sys.argv)
    feedback = f'Execution parameters:\t{para}'
    print(feedback, end="\n", flush=True)
    target_gene_table = read_table(complete_path(args.id_list))
    synteny_mapping_mutiprocess(target_gene_table, complete_path(args.bed), complete_path(args.out_dir),
                                args.size, complete_path(args.SynNet_f),
                                complete_path(args.species_namelist_file), args.threads)


if __name__ == "__main__":
    main()

