# CEG
# haochenli
# Python3
# v0.2.5
import os
import sys
import shutil
import subprocess
import re
import datetime
import numpy as np
import pandas as pd
import networkx as nx
import argparse
import igraph as ig
from collections import Counter
from sklearn.feature_extraction.text import CountVectorizer
from Bio import SeqIO


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


def check_external_software():
    # Ensure that the software in the environment is correctly configured
    required_software = ['exec_annotation']
    missing = [pkg for pkg in required_software if shutil.which(pkg) is None]
    if missing:
        raise EnvironmentError(f"Missing required software: {', '.join(missing)}")


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


def clean_network_id(clean_network):
    # Ensure that the correct ID is entered
    clean_network_table = read_table(clean_network)
    node_table = []
    for edge in clean_network_table:
        edge_sp = edge.split("\t")
        for node in edge_sp:
            node_table.append(node)
    node_set = set(node_table)
    return node_set


def get_block_id_from_edge(clean_network, net_f, retain_if):
    # Identify the IDs of the blocks to be extracted from the given edge file
    matched_block_id_table = []
    task_volume = enumerate_count(net_f)
    i = 0
    node_set = clean_network_id(clean_network)
    net_f_table = read_table(net_f)
    for l1 in net_f_table:
        i += 1
        print('find' + ' ' + str(i) + ' of ' + str(task_volume), end="\n", flush=True)
        l1_sp = l1.split("\t")
        l1_sp2 = l1_sp[0].split("-")
        if retain_if is True:
            if l1_sp[2] in node_set or l1_sp[3] in node_set:
                matched_block_id_table.append(f"{l1_sp[2]}&{l1_sp[3]}\t{l1_sp2[0]}")
        else:
            if l1_sp[2] in node_set and l1_sp[3] in node_set:
                matched_block_id_table.append(f"{l1_sp[2]}&{l1_sp[3]}\t{l1_sp2[0]}")
    return matched_block_id_table


def get_block_id_from_namelist(namelist_file, net_f, retain_if):
    # Identify the IDs of the blocks to be extracted from the given list of IDs
    matched_block_id_table = []
    task_volume = enumerate_count(net_f)
    i = 0
    namelist_table = read_table(namelist_file)
    node_set = set(namelist_table)
    net_f_table = read_table(net_f)
    for l1 in net_f_table:
        i += 1
        print('find' + ' ' + str(i) + ' of ' + str(task_volume), end="\n", flush=True)
        l1_sp = l1.split("\t")
        l1_sp2 = l1_sp[0].split("-")
        if retain_if is True:
            if l1_sp[2] in node_set or l1_sp[3] in node_set:
                matched_block_id_table.append(f"{l1_sp[2]}&{l1_sp[3]}\t{l1_sp2[0]}")
        else:
            if l1_sp[2] in node_set and l1_sp[3] in node_set:
                matched_block_id_table.append(f"{l1_sp[2]}&{l1_sp[3]}\t{l1_sp2[0]}")
    return matched_block_id_table


def filter_block_len(raw_result, block_size, output):
    # Trim the block lengths according to the requirements
    with open(output, 'a+', encoding='utf-8') as f:
        table = read_table(raw_result)
        for line in table:
            if not len(line.strip()) == 0:
                line_sp = line.split("\t")
                line_sp2 = line_sp[0].split("&")
                target_id = line_sp2[0]
                block_id = line_sp[1]
                table2 = []
                syn_info = "\t".join(line_sp[2:])
                syn_info_sp = syn_info.split(">")
                for line2 in syn_info_sp[:-1]:
                    table2.append(line2.rstrip("\n"))
            num = 0
            table3 = []
            for line3 in table2:
                line3_sp = line3.split(">")
                for line4 in line3_sp:
                    num = num + 1
                    if line4.__contains__(target_id):
                        table3.append(num)
                        continue
            try:
                min_num = int(table3[0]) - int(block_size)
            except IndexError:
                continue
            max_num = int(table3[0]) + int(block_size)
            num2 = 0
            table4 = []
            for line5 in table2:
                num2 = num2 + 1
                if min_num <= num2 <= max_num:
                    table4.append(line5)
            new_str = line_sp[0] + '\t' + block_id + "\t" + ">".join(table4) + '>'
            print(new_str, file=f)


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


def key2value(search_table, dictionary, output):
    # Obtain the value from the specified dictionary based on the key
    with open(output, 'a+', encoding='utf-8') as f:
        for i2 in search_table:
            i2_sp = i2.split("\t")
            value = dictionary.get(i2_sp[1])
            new_line = i2 + '\t' + value
            print(new_line.rstrip("\n"), file=f)


def transfer(raw_collinearity, result_path):
    # Transform the display format of the block
    intermediate_table = read_table(raw_collinearity)
    with open(result_path, 'a+', encoding='utf-8') as f1:
        for line3 in intermediate_table:
            line3_sp = line3.split("\t")
            syn_info = "\t".join(line3_sp[2:])
            for line4 in syn_info:
                f1.write(line4.replace('>', '\n'))


def table2file(inputtable, outputfile):
    # Format the output data
    with open(outputfile, 'a+', encoding='utf-8') as f:
        for line in inputtable:
            print(line.rstrip("\n"), file=f)


def table2file2(inputtable, outputfile):
    # Format the output data
    with open(outputfile, 'a+', encoding='utf-8') as f:
        print('source\ttarget', file=f)
        for line in inputtable:
            str_line = str(line)
            line_sp = str_line.split("'")
            s1 = line_sp[1]
            s2 = line_sp[3]
            out_line = s1 + "\t" + s2
            print(out_line, file=f)


def table2file3(inputtable, outputfile):
    # Format the output data
    with open(outputfile, 'a+', encoding='utf-8') as f:
        print('nodes\tdegree', file=f)
        for line in inputtable:
            str_line = str(line)
            line_sp = str_line.split("'")
            s1 = line_sp[1]
            line_sp2 = str_line.split(" ")
            s2 = line_sp2[1].rstrip(")")
            out_line = s1 + "\t" + s2
            print(out_line, file=f)


def frozenset_progressive(filter_table, out_path1, out_path2):
    # Classify the nodes in the network using the greedy algorithm
    sub_syn_network = nx.Graph()
    for line in filter_table:
        str_line = str(line)
        str_line_sp = str_line.split("'")
        str1 = str_line_sp[1]
        str2 = str_line_sp[3]
        new_line = str1 + '\t' + str2
        head, tail = [str(x) for x in new_line.split()]
        sub_syn_network.add_edge(head, tail)
    frozenset_filter_table = nx.algorithms.community.greedy_modularity_communities(sub_syn_network)
    with open(out_path1, 'a+', encoding='utf-8') as f3:
        with open(out_path2, 'a+', encoding='utf-8') as f4:
            num = 0
            for line4 in frozenset_filter_table:
                line4_str = str(line4)
                if line4_str.startswith("frozenset"):
                    num = num + 1
                    line4_sp = line4_str.split("{")
                    line4_sp2 = line4_sp[1].split("}")
                    line4_sp3 = line4_sp2[0]
                    line4_rp1 = line4_sp3.replace("'", "").replace(" ", "").replace(",", "\t")
                    new_line4 = str(num) + '\t' + line4_rp1
                    print(new_line4.rstrip("\n"), file=f3)
                    line4_sp4 = line4_rp1.split("\t")
                    for line5 in line4_sp4:
                        new_line5 = line5 + '\t' + str(num)
                        print(new_line5.rstrip("\n"), file=f4)


def k_core_progressive(alignment_file, output_path, base_filename):
    # Filter the core network based on the k-core algorithm
    syn_network = nx.Graph()
    with open(alignment_file) as file:
        for line in file:
            head, tail = [str(x) for x in line.split()]
            syn_network.add_edge(head, tail)

        content = 1
        num = 0
        while content == 1:
            num += 1
            out_path1 = f"{output_path}/nodes/{base_filename}_k{str(num)}.nodes"
            out_path2 = f"{output_path}/edges/{base_filename}_k{str(num)}.edges.tsv"
            out_path3 = f"{output_path}/degree/{base_filename}_k{str(num)}.degree"
            out_path4 = f"{output_path}/group/{base_filename}_k{str(num)}.group"
            out_path5 = f"{output_path}/group/{base_filename}_k{str(num)}.group.2col"
            filter_table1 = nx.k_core(syn_network, num).nodes
            if len(filter_table1) == 0:
                content = 0
            else:
                print('calculating k-core k=' + str(num), end="\n", flush=True)
                table2file(filter_table1, out_path1)
                print('calculating edges  k=' + str(num), end="\n", flush=True)
                filter_table2 = nx.k_core(syn_network, num).edges
                table2file2(filter_table2, out_path2)
                print('calculating degree k=' + str(num), end="\n", flush=True)
                filter_table3 = nx.k_core(syn_network, num).degree
                table2file3(filter_table3, out_path3)
                print('calculating group k=' + str(num), end="\n", flush=True)
                frozenset_progressive(filter_table2, out_path4, out_path5)


def standard_block_output(input_file, standard_block_path):
    # Output a file with block information in standard format
    with open(standard_block_path, 'a+', encoding='utf-8') as f:
        table1 = read_table(input_file)
        table2 = []
        for line2 in table1:
            line2_sp = line2.split("\t")
            text = 'Target_gene:>' + line2_sp[0] + '>' + 'Block_id:>' + line2_sp[1] + \
                   ">Synteny_info:>" + "\t".join(line2_sp[2:])
            table2.append(text)
        table2.sort()
        table2 = list(set(table2))
        for line2 in table2:
            new_text = line2.replace('>', '\n')
            print(new_text, file=f)


def extract_new(fasta_file, wanted_file, result_file):
    # Extract protein sequences based on IDs.
    wanted = set()
    with open(wanted_file) as f:
        for line in f:
            line = line.strip()
            if line != "":
                wanted.add(line)

    fasta_sequences = SeqIO.parse(open(fasta_file), 'fasta')
    with open(result_file, "w") as f:
        for seq in fasta_sequences:
            if seq.id in wanted:
                SeqIO.write([seq], f, "fasta")


def build_graph(edge_file_input):
    # Create a network with networkX
    g1 = nx.Graph()
    edge_table = read_table(edge_file_input)
    for l1 in edge_table:
        l1_sp = l1.split("\t")
        g1.add_edge(l1_sp[0], l1_sp[1])
    return g1


def get_id_dic(edge_file_input):
    # Obtain node from an edge file
    prefix_table = read_table(edge_file_input)
    id_table = []
    id_dic = {}
    for l1 in prefix_table:
        l1_sp = l1.split("\t")
        for l2 in l1_sp:
            id_table.append(l2)
    prefix_list2 = sorted(list(set(id_table)))
    count1 = -1
    for l3 in prefix_list2:
        count1 = count1 + 1
        id_dic[count1] = l3
    return id_dic


def new_id_file_generation(edge_file_input, id_dic, new_edge_file_input):
    # Translate text IDs into pure numeric IDs
    id_dic_reverse = {v: k for k, v in id_dic.items()}
    with open(new_edge_file_input, 'a+', encoding='utf-8') as f1:
        edge_table = read_table(edge_file_input)
        for l4 in edge_table:
            if not l4.startswith("source"):
                l4_sp = l4.split("\t")
                str1 = l4.replace(l4_sp[0], str(id_dic_reverse.get(l4_sp[0])))
                str2 = str1.replace(l4_sp[1], str(id_dic_reverse.get(l4_sp[1])))
                print(str2, file=f1)


def classify_community_infomap(new_edge_file_input):
    # Classify the nodes in the network using the Infomap algorithm
    df = pd.read_table(new_edge_file_input, sep='\t', names='source\ttarget'.split('\t'))
    graph1 = ig.Graph.DataFrame(df)
    out_table = ig.Graph.community_infomap(graph1)
    return out_table


def get_infomap_col2(new_edge_file, id_dic, infomap_classify_col2_file):
    # Translate the Infomap clustering results into a two-column file
    infomap_tuple = classify_community_infomap(new_edge_file)
    with open(infomap_classify_col2_file, 'a+', encoding='utf-8') as f1:
        count2 = -1
        for l5 in infomap_tuple:
            count2 = count2 + 1
            for l6 in l5:
                out_str1 = id_dic.get(l6) + '\t' + str(count2)
                print(out_str1, file=f1)


def get_size_frequency(new_edge_file, size_frequency_input):
    # Calculate the size and frequency of each group under Infomap classification
    infomap_tuple = classify_community_infomap(new_edge_file)
    with open(size_frequency_input, 'a+', encoding='utf-8') as f1:
        print("size" + "\t" + "count" + "\t" + "frequency", file=f1)
        size_table = []
        counted_size = {}
        count3 = -1
        for l7 in infomap_tuple:
            count3 = count3 + 1
            size_table.append(len(l7))
        for size in size_table:
            counted_size[size] = counted_size.get(size, 0) + 1
        count_sum = 0
        for count in counted_size.values():
            count_sum = count_sum + count
        for l8 in counted_size.keys():
            size = counted_size.get(l8)
            fq = size / count_sum
            out_str2 = str(l8) + "\t" + str(size) + "\t" + str(fq)
            print(out_str2, file=f1)


def infomap_clustering(edge_file, new_edge_file, infomap_classify_col2, size_frequency_file):
    # Implementation of the Infomap clustering function
    id_dic = get_id_dic(edge_file)
    new_id_file_generation(edge_file, id_dic, new_edge_file)
    get_infomap_col2(new_edge_file, id_dic, infomap_classify_col2)
    get_size_frequency(new_edge_file, size_frequency_file)


def match_species_classification_information(k1_node, sp_info_file, group_2col, infomap_2col, output_file):
    # Match IDs with their corresponding various classification information
    sp_info_table = read_table(sp_info_file)
    node_table = read_table(k1_node)
    cls_info_table = []
    with open(output_file, 'a+', encoding='utf-8') as f1:
        print("#ID\tClade\tOrder\tFamily\tSpecies_Name\tGreedy_cluster\tInfomap_cluster", file=f1)
        for l1 in sp_info_table:
            if not l1.startswith("#"):
                l1_sp = l1.split("\t")
                cls_info_table.append(l1_sp[1:])
        sp_info_dic = {}
        for info in cls_info_table:
            key = info[0]
            value = "\t".join(info[1:])
            sp_info_dic[key] = value
        group_dic = read_dic(group_2col)
        infomap_dic = read_dic(infomap_2col)
        for l2 in node_table:
            cls = sp_info_dic.get(l2[0:4]).strip("\n")
            group_id = group_dic.get(l2).strip("\n")
            infomap_id = infomap_dic.get(l2).strip("\n")
            new_str = f"{l2}\t{cls}\t{group_id}\t{infomap_id}"
            print(new_str, file=f1)


def get_file_uniq_title(file_name):
    # Backup, to prevent IDs with special prefixes from occurring
    file_name_sp = file_name.split('/')
    title = file_name_sp[-1]
    if title.startswith('Alp'):
        uniq_title = title[0:3]
    elif title.startswith('XS'):
        uniq_title = title[0:3]
    else:
        uniq_title = title[0:4]
    return uniq_title


def create_prefix_file(namelist):
    # Backup，to prevent IDs with special prefixes from occurring
    all_loc_talbe = []
    file_list = read_table(namelist)
    for l1 in file_list:
        file_table = read_table(l1)
        new_file_name = l1 + '.new'
        with open(new_file_name, 'a+', encoding='utf-8') as f:
            for l2 in file_table:
                new_line = get_file_uniq_title(l1) + '\t' + l2
                all_loc_talbe.append(new_line)
                print(new_line.rstrip("\t"), file=f)
    return all_loc_talbe


def get_data_frame_index(namelist):
    # Retrieve the index
    file_list = read_table(namelist)
    table = []
    for l1 in file_list:
        if not l1.startswith("#"):
            l1_sp = l1.split("\t")
            table.append(l1_sp[5])
    return table


def get_loc_table(count_table, sp_cls_table):
    # Retrieve the position and value
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


def block_num_stat(block_file):
    # Count the number of blocks generated by pairwise comparisons between species in the dataset that make up the SNN
    output_table = []
    stat_table = []
    block_table = read_table(block_file)
    block_set = set()
    for l1 in block_table:
        l1_sp = l1.split("\t")
        block_set.add(l1_sp[1])
    for l2 in block_set:
        text = re.sub(r'\d+', '', l2)
        stat_table.append(text)
    count = Counter(stat_table)
    for element, frequency in count.items():
        output_table.append(f"{element}\t{frequency}")
    return output_table


def block_length_sum_stat(block_file):
    """
    Calculate the total length of blocks generated by pairwise comparisons between species in
    the dataset that constitute the SNN
    """
    output_table = []
    dic_table = []
    total_counts = {}
    block_table = read_table(block_file)
    block_set = set()
    for l1 in block_table:
        l1_sp = l1.split("\t")
        block_set.add("\t".join(l1_sp[1:]))
    for l2 in block_set:
        l2_sp = l2.split("\t")
        text = re.sub(r'\d+', '', l2_sp[0])
        block_info = "\t".join(l2_sp[1:])
        flag = ">"
        counter = Counter(block_info)
        count = counter[flag]
        tem_table = [text, count]
        dic_table.append(tem_table)
    for entry, count in dic_table:
        if entry in total_counts:
            total_counts[entry] += count
        else:
            total_counts[entry] = count
    for entry, total in total_counts.items():
        output_table.append(f"{entry}\t{total}")
    return output_table


def merge_pep(input_namelist, data_path, output_p):
    # Merge the PEP files of species that constitute the SNN
    species_table_1col = []
    species_table = read_table(input_namelist)
    for specie in species_table:
        if not specie.startswith("#"):
            specie_sp = specie.split("\t")
            species_table_1col.append(specie_sp[0])
    for species in species_table_1col:
        with open(f"{output_p}/all.pep", 'a+', encoding='utf-8') as f1:
            pep_table = read_table(f"{data_path}/{species}.pep")
            for l1 in pep_table:
                print(l1.rstrip("\n"), file=f1)


def fill_dataframe_num(loc_info, dataframe):
    # Fill in the matrix based on the statistical results
    for l1 in loc_info:
        if len(l1) == 2:
            row_loc = l1[0]
            col_loc = l1[0]
            value = l1[1]
            dataframe.loc[row_loc, col_loc] = int(value)
        if len(l1) == 3:
            row_loc = l1[0]
            col_loc = l1[1]
            value = l1[2]
            dataframe.loc[row_loc, col_loc] = int(value)
            dataframe.loc[col_loc, row_loc] = int(value)


def block_stat(namelist_file_input, filtered_block_file_input, out_dir):
    # Implementation of the --block_stat function
    os.mkdir(f"{out_dir}/block_stat")
    data_frame_col = get_data_frame_index(namelist_file_input)
    data_frame_row = get_data_frame_index(namelist_file_input)
    block_num_stat_table = block_num_stat(filtered_block_file_input)
    block_length_sum_stat_table = block_length_sum_stat(filtered_block_file_input)
    df_num = pd.DataFrame(columns=data_frame_col, index=data_frame_row)
    pd.set_option('future.no_silent_downcasting', True)
    df_num = df_num.fillna(0).infer_objects(copy=False)
    df_len_sum = pd.DataFrame(columns=data_frame_col, index=data_frame_row)
    df_len_sum = df_len_sum.fillna(0).infer_objects(copy=False)
    fill_dataframe_num(get_loc_table(block_num_stat_table, namelist_file_input), df_num)
    fill_dataframe_num(get_loc_table(block_length_sum_stat_table, namelist_file_input), df_len_sum)
    df_len_avg = df_len_sum / df_num
    df_num.to_csv(f"{out_dir}/block_stat/block_num_stat.tsv",
                  sep='\t', index=True, header=True)
    df_len_sum.to_csv(f"{out_dir}/block_stat/block_len_sum_stat.tsv",
                      sep='\t', index=True, header=True)
    df_len_avg.to_csv(f"{out_dir}/block_stat/block_len_avg_stat.tsv",
                      sep='\t', index=True, header=True)


def add_kegg_info(kofamscan_output, out_path1, sp_cluster_file, out_path2):
    # Merge KEGG annotation results with classification information
    kofamscan_table = read_table(kofamscan_output)
    with open(out_path1, 'a+', encoding='utf-8') as f1:
        for l1 in kofamscan_table:
            if l1.startswith("*"):
                l1_sp = l1.split("\t")
                gene_name = l1_sp[1]
                ko_id = l1_sp[2]
                ko_definition = l1_sp[6].replace("\"", "")
                new_text = f"{gene_name}\t{ko_id}\t{ko_definition}"
                print(new_text.rstrip("\n"), file=f1)
    ko_dic = read_dic(out_path1)
    sp_cluster_table = read_table(sp_cluster_file)
    with open(out_path2, 'a+', encoding='utf-8') as f2:
        print("id\tclade\torder\tfamily\tspecies_name\tKO_id\tKO_definition\tgreedy_cluster\tinfomap_cluster", file=f2)
        for l2 in sp_cluster_table:
            if not l2.startswith("#"):
                l2_sp = l2.split("\t")
                try:
                    match_value = ko_dic.get(l2_sp[0]).rstrip("\n")
                except AttributeError:
                    match_value = "None\tNone"
                out_str1 = "\t".join(l2_sp[:5])
                out_str2 = "\t".join(l2_sp[-2:])
                new_text = f"{out_str1}\t{match_value}\t{out_str2}"
                print(new_text.rstrip("\n"), file=f2)


def block_info_cluster_info_trans(input_cluster_col2, input_syn_clean_result, out_put):
    # Convert gene IDs in the block into classification information
    infomap_dic = read_dic(input_cluster_col2)
    with open(out_put, 'a+', encoding='utf-8') as f:
        clean_result_table = read_table(input_syn_clean_result)
        for l1 in clean_result_table:
            l1_sp = l1.split("\t")
            line_table = []
            for l2 in l1_sp[2:]:
                l2_sp = l2.split(">")
                for l3 in l2_sp:
                    if infomap_dic.get(l3) is not None:
                        line_table.append(int(infomap_dic.get(l3)))
            line_table = list(set(line_table))
            out_str = '\t'.join(map(str, line_table))
            print(f"{l1_sp[0]}\t{l1_sp[1]}\t{out_str}", file=f)


def block_info_countvectorizer(input_block_info_cluster_file, block_info_df_out):
    # Convert the classification information of the block into a feature matrix
    block_info_cluster_table = read_table(input_block_info_cluster_file)
    block_info_cluster_array = []
    id_block_table = []
    vocab = []
    for l1 in block_info_cluster_table:
        l1_sp = l1.split("\t")
        line_table = [l1_sp[0], l1_sp[1]]
        id_block_table.append(line_table)
        block_info_cluster_array.append("\t".join(l1_sp[2:]))
        for l2 in l1_sp[2:]:
            vocab.append(l2)
    df_id = pd.DataFrame(id_block_table, columns=['id', 'block_id'])
    vocab = (set(vocab))
    vectorizer = CountVectorizer(vocabulary=vocab)
    count_matrix = vectorizer.fit_transform(block_info_cluster_array)
    count_array = count_matrix.toarray()
    df_info = pd.DataFrame(count_array)
    df_all = df_id.join(df_info)
    df_all.to_csv(block_info_df_out, sep="\t", index=False, header=False)


def block_specificity_score(block_info_df, out_put):
    # Calculate the specificity score for each block
    matrix_table = read_table(block_info_df)
    matrix_prefix_table = []
    block_id_table = []
    target_id_table = []
    for l1 in matrix_table:
        l1_sp = l1.split("\t")
        target_id_table.append(l1_sp[0])
        block_id_table.append(l1_sp[1])
        line_table = []
        for l2 in l1_sp[2:]:
            line_table.append(l2)
        matrix_prefix_table.append(line_table)

    matrix = np.array(matrix_prefix_table)
    arr_numeric = matrix.astype(float)
    # 计算每一列的平均值
    col_means = np.mean(arr_numeric, axis=0)
    # 计算每一列的标准差
    col_stds = np.std(arr_numeric, axis=0)
    col_stds[col_stds == 0] = 1e-8
    # 标准化处理，然后取绝对值
    standardized_matrix = np.abs((arr_numeric - col_means) / col_stds)
    # 计算每一行的和
    row_sums = np.sum(standardized_matrix, axis=1)
    count = 0
    with open(out_put, 'a+', encoding='utf-8') as f1:
        names = 'target_gene\tblock_id\tblock_specificity_score'
        print(names.rstrip("\n"), file=f1)
        for l3 in range((len(block_id_table)) - 1):
            count = count + 1
            out_str = f"{target_id_table[count]}\t{block_id_table[count]}\t{row_sums[count]}"
            print(out_str.rstrip("\n"), file=f1)


def main():
    text1 = '''The script will collect the upstream and downstream genes in the collinearity blocks for 
        the genes listed.; 
           '''

    text2 = '''-l: list of gene ID; input a list of gene ID;
        The option -I can be selected to input an edge file as an alternative.
        '''

    text3 = '''-e: edge_file; input a edge file with two columns;
        '''

    text4 = '''-b: input SynNet_f file; 
        Input SynNet_f file which generated through synetfind.py script
        '''

    text5 = '''-N: input SynNet.prefix file; 
        Input SynNet.prefix file which generated through synetprefix.py script
        '''

    text6 = '''-S: determine the block size (default=10); Determine the maximum number of genes to be retained upstream 
        and downstream of the target gene in the collinearity block.
        '''

    text7 = '''-d: directory of pep_file and bed_file; 
       Directory where the pep file and bed file are located
       '''

    text8 = '''-i: list of species names; 
        input a list of file names for each species' pep and bed files
        '''

    text9 = '''-p: number of threads: Number of threads for functional annotation (default=8)
        '''

    text10 = '''-r: retain additional results; 
        Retain the matching results of non listed genes
        '''

    text11 = '''-o: output_dir; 
        Output files to the specified directory
        '''

    text12 = '''--KEGG: KEGG annotation; 
        Annotate nodes in a network using kofamsacn
        '''

    text13 = '''--block_stat: Statistically analyze the collinear blocks that form the network
        '''

    text14 = '''-P: input all species pep_file 
        '''

    text15 = '''--block_mining: Create a block feature matrix and calculate the specificity score 
        '''

    parser = argparse.ArgumentParser(description=text1)
    group = parser.add_mutually_exclusive_group(required=True)
    group2 = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument('-i', '--species_namelist_file', type=str, required=True, help=text8)
    group.add_argument('-l', '--id_list', type=str, help=text2)
    group.add_argument('-e', '--edge_file', type=str, help=text3)
    parser.add_argument('-b', '--SynNet_f', type=str, required=True, help=text4)
    parser.add_argument('-N', '--SynNet_prefix', type=str, required=True, help=text5)
    group2.add_argument('-d', '--data_path', type=str, help=text7)
    group2.add_argument('-P', '--all_pep', type=str, help=text14)
    parser.add_argument('-o', '--out_dir', type=str, required=True, help=text11)
    parser.add_argument('-S', '--size', type=int, required=True, default=10, help=text6)
    parser.add_argument('-p', '--threads', type=int, default=8, help=text9)
    parser.add_argument('-r', '--retain', action='store_true', help=text10)
    parser.add_argument('--KEGG', action='store_true', help=text12)
    parser.add_argument('--block_stat', action='store_true', help=text13)
    parser.add_argument('--block_mining', action='store_true', help=text15)
    args = parser.parse_args()
    para = " ".join(sys.argv)
    feedback = f'Execution parameters:\t{para}'
    print(feedback, end="\n", flush=True)

    block_id_table = []
    base_filename = ''
    foldername = "SNN_" + datetime.datetime.now().strftime("%Y%m%d_%H%M")

    output_path = f"{complete_path(args.out_dir)}/{foldername}"
    os.mkdir(output_path)
    path_sep = "/"
    if args.retain:
        retain = True
    else:
        retain = False
    if args.id_list:
        base_filename = os.path.basename(args.id_list)
        block_id_table = get_block_id_from_namelist(complete_path(args.id_list), complete_path(args.SynNet_f), retain)
    if args.edge_file:
        base_filename = os.path.basename(args.edge_file)
        block_id_table = get_block_id_from_edge(complete_path(args.edge_file), complete_path(args.SynNet_f), retain)
    raw_result_path = output_path + path_sep + 'raw_result'
    result_path = output_path + path_sep + 'result'
    os.makedirs(output_path + path_sep + 'K_core_result')
    os.makedirs(output_path + path_sep + 'K_core_result' + path_sep + 'nodes')
    os.makedirs(output_path + path_sep + 'K_core_result' + path_sep + 'edges')
    os.makedirs(output_path + path_sep + 'K_core_result' + path_sep + 'degree')
    os.makedirs(output_path + path_sep + 'K_core_result' + path_sep + 'group')
    k_core_result_path = output_path + path_sep + 'K_core_result'
    raw_result_filter_path = output_path + path_sep + 'raw_result_filter'
    k1_node_path = k_core_result_path + path_sep + 'nodes' + path_sep + base_filename + '_''k' + '1' + '.nodes'
    k1_edge_path = k_core_result_path + path_sep + 'edges' + path_sep + base_filename + '_''k' + '1' + '.edges.tsv'
    new_k1_edge_file = k_core_result_path + path_sep + 'edges' + path_sep + base_filename + '_''k' + '1' + '.num_ID.tsv'
    k1_group_path = k_core_result_path + path_sep + 'group' + path_sep + base_filename + '_''k' + '1' + '.group.2col'
    infomap_classify_2col_file = output_path + path_sep + base_filename + '_''k' + '1' + '_infomap.2col.tsv'
    size_frequency_file = output_path + path_sep + base_filename + '_''k' + '1' + '_infomap.size_fq.tsv'
    sp_cluster_sum_file = output_path + path_sep + base_filename + '_''k' + '1' + '_sp_cluster.tsv'
    sp_ko_cls_sum_file = output_path + path_sep + base_filename + '_''k' + '1' + '_sp_KO_cluster.tsv'
    standard_block_path = output_path + path_sep + 'block_info'
    print('reading dic', end="\n", flush=True)
    synnet_prefix_dic = read_dic(complete_path(args.SynNet_prefix))
    if args.data_path:
        merge_pep(args.species_namelist_file, complete_path(args.data_path), output_path)
    print('Exporting results')
    key2value(block_id_table, synnet_prefix_dic, raw_result_path)
    print('Block size filtering', end="\n", flush=True)
    filter_block_len(raw_result_path, args.size, raw_result_filter_path)
    print('Expanding files', end="\n", flush=True)
    transfer(raw_result_filter_path, result_path)
    print('Generating list', end="\n", flush=True)
    standard_block_output(raw_result_filter_path, standard_block_path)
    print('Calculating network', end="\n", flush=True)
    k_core_progressive(result_path, k_core_result_path, base_filename)
    print('Infomap clustering', end="\n", flush=True)
    infomap_clustering(k1_edge_path, new_k1_edge_file, infomap_classify_2col_file, size_frequency_file)
    print('Match species classification information', end="\n", flush=True)
    match_species_classification_information(k1_node_path, complete_path(args.species_namelist_file), k1_group_path,
                                             infomap_classify_2col_file, sp_cluster_sum_file)
    print('Extract protein sequences from the network', end="\n", flush=True)
    if args.data_path:
        extract_new(f"{output_path}/all.pep", k1_node_path, f"{output_path}/{os.path.basename(k1_node_path)}.pep")
    if args.all_pep:
        extract_new(args.all_pep, k1_node_path, f"{output_path}/{os.path.basename(k1_node_path)}.pep")
    if args.block_stat:
        print('Calculate block information', end="\n", flush=True)
        block_stat(complete_path(args.species_namelist_file), raw_result_filter_path, output_path)
    if args.block_mining:
        print('Calculate block specificity score', end="\n", flush=True)
        os.makedirs(output_path + path_sep + 'specific_block')
        block_info_cluster_file = f'{output_path}/specific_block/{base_filename}_k1_block_info.tsv'
        block_info_df = f'{output_path}/specific_block/{base_filename}_k1_block_info_matrix.tsv'
        block_info_score = f'{output_path}/specific_block/{base_filename}_k1_block_info_score.tsv'
        block_info_cluster_info_trans(k1_group_path, raw_result_filter_path, block_info_cluster_file)
        block_info_countvectorizer(block_info_cluster_file, block_info_df)
        block_specificity_score(block_info_df, block_info_score)

    if args.KEGG:
        check_external_software()
        print('Start KEGG annotation process')
        os.mkdir(f"{output_path}/KEGG_annotation")
        subprocess.run(["exec_annotation", "-o",
                        f"{output_path}/KEGG_annotation/{os.path.basename(k1_node_path)}.KO.tsv",
                        f"{output_path}/{os.path.basename(k1_node_path)}.pep",
                        f"--cpu={args.threads}", "-f", "detail-tsv", f"--tmp-dir={output_path}"], check=True)
        add_kegg_info(f"{output_path}/KEGG_annotation/{os.path.basename(k1_node_path)}.KO.tsv",
                      f"{output_path}/KEGG_annotation/{os.path.basename(k1_node_path)}.KO_k1.tsv",
                      sp_cluster_sum_file, sp_ko_cls_sum_file)


if __name__ == "__main__":
    main()
