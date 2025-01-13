# CEG
# haochenli
# Python3
# v0.2.5
import glob
import shutil
import subprocess
import os
import sys
import datetime
import argparse
import re
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


def complete_path(input_path):
    # Ensure that the input is an absolute path
    if not os.path.isabs(input_path):
        input_path = os.path.join(os.getcwd(), input_path)
    return input_path


def check_external_software():
    # Ensure that the software in the environment is correctly configured
    required_software = ['hmmsearch']
    missing = [pkg for pkg in required_software if shutil.which(pkg) is None]
    if missing:
        raise EnvironmentError(f"Missing required software: {', '.join(missing)}")


def extract_new(fasta_file, wanted_file, result_file):
    # Extract the data
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


def merge_files(output_filename, input_filenames):
    # Merge files
    with open(output_filename, 'w') as outfile:
        for filename in input_filenames:
            with open(filename, 'r') as infile:
                outfile.write(infile.read())


def find_record(input_namelist, input_hmm, data_path, output_dir, evalue, threads, synet_file, retain_if):
    # Search for sequences based on the given HMM file and extract the files where the HMM model hits
    check_external_software()
    input_hmm_name = os.path.basename(input_hmm)
    # Define the output directory
    foldername = "SynNet_" + input_hmm_name + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M")
    result_path = f"{output_dir}/{foldername}"
    os.makedirs(result_path)
    os.chdir(result_path)
    species_table_1col = []
    species_table = read_table(input_namelist)
    for specie in species_table:
        if not specie.startswith("#"):
            specie_sp = specie.split("\t")
            species_table_1col.append(specie_sp[0])
    # Perform sequence retrieval using hmmsearch
    for species in species_table_1col:
        print(f"Searching proteins by module:{input_hmm_name}   sp_abb:{species}", end="\n", flush=True)
        with open(f"{result_path}/{species}.hmmout", 'a+', encoding='utf-8') as f1:
            subprocess.run(["hmmsearch", "--noali", "-E", str(evalue), "--cpu", str(threads),
                            f"{input_hmm}", f"{data_path}/{species}.pep"],
                           check=True, stdout=f1)
            hmmout_table = read_table(f"{result_path}/{species}.hmmout")
        with open(f"{result_path}/{species}.hmm_genelist_f", 'a+', encoding='utf-8') as f2:
            print(f"1iE-value\tscore\tbias\tE-value\tscore\tbias\texp\tN\tSequence", file=f2)
            for l1 in hmmout_table[19:]:
                if l1.strip() != '':
                    new_str = re.sub(r'\s+', ' ', l1).rstrip(" ")
                    str_sp1 = new_str.split(" ")
                    out_str = "\t".join(str_sp1)
                    print(out_str, file=f2)
                else:
                    break

        # Organize the results of hmmsearch to determine the species IDs where the HMM model has hits
        with open(f"{result_path}/{species}.hmm_genelist", 'a+', encoding='utf-8') as f3:
            for l2 in hmmout_table:
                if l2.startswith(">>"):
                    l2_sp = l2.split(" ")
                    print(l2_sp[1], file=f3)
        extract_new(f"{data_path}/{species}.pep", f"{result_path}/{species}.hmm_genelist",
                    f"{result_path}/{species}.hmm_genelist.faa")
    print(f"merge genelist", end="\n", flush=True)
    merge_files(f"{result_path}/{input_hmm_name}.genes", glob.glob(f"{result_path}/*.hmm_genelist"))
    merge_files_set = set(read_table(f"{result_path}/{input_hmm_name}.genes"))

    # Combine the IDs of hmmsearch hits from various species
    with open(f"{result_path}/all.hmm_genelist", 'a+', encoding='utf-8') as f3_2:
        for l2_2 in merge_files_set:
            print(l2_2, file=f3_2)

    # Extract corresponding results from the overall synteny network based on the list
    synet_table = read_table(synet_file)
    synet_table_match = []
    with open(f"{result_path}/{input_hmm_name}.genelist_SynNet_f", 'a+', encoding='utf-8') as f4:
        for l3 in synet_table:
            l3_set = l3.strip().split("\t")
            if merge_files_set.intersection(l3_set):
                synet_table_match.append(l3)
                print(l3, file=f4)

    if retain_if is True:
        with open(f"{result_path}/{input_hmm_name}.genelist_SynNet_f_2col", 'a+', encoding='utf-8') as f5:
            for l4 in synet_table_match:
                l4_sp = l4.split("\t")
                out_str = f"{l4_sp[2]}\t{l4_sp[3]}"
                print(out_str, file=f5)

    else:
        with open(f"{result_path}/{input_hmm_name}.genes_cleaned-network", 'a+', encoding='utf-8') as f6:
            for l4 in synet_table_match:
                l4_sp = l4.split("\t")
                if l4_sp[2] in merge_files_set and l4_sp[3] in merge_files_set:
                    out_str = f"{l4_sp[2]}\t{l4_sp[3]}"
                    print(out_str, file=f6)


def find_record_customized(input_customized_namelist, output_dir, synet_file, retain_if):
    # Extract information from the overall synteny network based on a customized list
    check_external_software()
    customized_name = os.path.basename(input_customized_namelist)
    foldername = "SynNet_" + customized_name + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M")
    result_path = f"{output_dir}/{foldername}"
    os.makedirs(result_path)
    os.chdir(result_path)

    synet_table = read_table(synet_file)
    synet_table_match = []
    merge_files_set = set(read_table(input_customized_namelist))
    with open(f"{result_path}/{customized_name}.genelist_SynNet_f", 'a+', encoding='utf-8') as f4:
        for l3 in synet_table:
            l3_set = l3.strip().split("\t")
            if merge_files_set.intersection(l3_set):
                synet_table_match.append(l3)
                print(l3, file=f4)

    if retain_if is True:
        with open(f"{result_path}/{customized_name}.genelist_SynNet_f_2col", 'a+', encoding='utf-8') as f5:
            for l4 in synet_table_match:
                l4_sp = l4.split("\t")
                out_str = f"{l4_sp[2]}\t{l4_sp[3]}"
                print(out_str, file=f5)

    else:
        with open(f"{result_path}/{customized_name}.genes_cleaned-network", 'a+', encoding='utf-8') as f6:
            for l4 in synet_table_match:
                l4_sp = l4.split("\t")
                if l4_sp[2] in merge_files_set and l4_sp[3] in merge_files_set:
                    out_str = f"{l4_sp[2]}\t{l4_sp[3]}"
                    print(out_str, file=f6)


def main():
    global retain
    text1 = '''The script searches for corresponding nodes in the complete 
        network and constructs a subnetwork based on the HMM model.; 
           Usage: python synetfind.py -i species_namelist_file -d data_path -m module.hmm -n SynNet_file -o output_ptah
           '''
    text2 = '''-i: list of species names; 
        input a list of file names for each species' pep and bed files
        '''
    text3 = '''-m: input hmm_file
        '''
    text4 = '''-l: input customized id list
        '''
    text5 = '''-d: directory of pep_file and bed_file; 
       Directory where the pep file and bed file are located
       '''
    text6 = '''-n: input Synnet_file; 
        input Synnet_file which generated through SynBuild
        '''
    text7 = '''-o: output_dir; 
       Output files to the specified directory
       '''
    text8 = '''-E: hmmsearch e-value; (default=0.001)
        '''
    text9 = '''-p: number of threads: Used for hmmsearch (default=8)'''
    text10 = '''-r: retain additional results; 
        IDs that have not been hit in hmmssearch will also be retained
        '''

    parser = argparse.ArgumentParser(description=text1)
    group = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument('-i', '--species_namelist_file', type=str, required=False, help=text2)
    group.add_argument('-m', '--hmm_file', type=str, help=text3)
    group.add_argument('-l', '--custom_id_list', type=str, help=text4)
    parser.add_argument('-d', '--data_path', type=str, required=False, help=text5)
    parser.add_argument('-n', '--synnet_file', type=str, required=True, help=text6)
    parser.add_argument('-o', '--out_dir', type=str, required=True, help=text7)
    parser.add_argument('-E', '--evalue', type=float, default=0.001, help=text8)
    parser.add_argument('-p', '--threads', type=int, default=8, help=text9)
    parser.add_argument('-r', '--retain', action='store_true', help=text10)
    args = parser.parse_args()
    para = " ".join(sys.argv)
    feedback = f'Execution parameters:\t{para}'
    print(feedback, end="\n", flush=True)
    if args.retain:
        retain = True
    else:
        retain = False
    if args.hmm_file:
        find_record(complete_path(args.species_namelist_file), complete_path(args.hmm_file),
                    complete_path(args.data_path),
                    complete_path(args.out_dir), args.evalue, args.threads,
                    complete_path(args.synnet_file), retain)
    if args.custom_id_list:
        find_record_customized(complete_path(args.custom_id_list), complete_path(args.out_dir),
                               complete_path(args.synnet_file), retain)


if __name__ == "__main__":
    main()
