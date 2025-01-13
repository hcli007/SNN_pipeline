# CEG
# haochenli
# Python3
# v0.2.5
import subprocess
import os
import sys
import datetime
import time
import shutil
import itertools
import argparse


def check_external_software():
    # Ensure that the software in the environment is correctly configured
    required_software = ['MCScanX', 'diamond']
    missing = [pkg for pkg in required_software if shutil.which(pkg) is None]
    if missing:
        raise EnvironmentError(f"Missing required software: {', '.join(missing)}")


def complete_path(input_path):
    # Ensure that the input is an absolute path
    if not os.path.isabs(input_path):
        input_path = os.path.join(os.getcwd(), input_path)
    return input_path


def merge_files(output_filename, input_filenames):
    # Merge files
    with open(output_filename, 'w') as outfile:
        for filename in input_filenames:
            with open(filename, 'r') as infile:
                outfile.write(infile.read())


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


def synet_build(input_file, data_path, output_path, hits, anchors, gaps, threads, duplicate, tandem):
    # Construct the overall synteny network for the species
    check_external_software()

    # Ensure that the server currently has sufficient computational resources
    while True:
        try:
            running = subprocess.check_output(["jobs", "-p"], universal_newlines=True)
            if len(running.strip().split('\n')) < int(threads):
                break
            time.sleep(1)
        except FileNotFoundError:
            pass

    species_table = read_table(input_file)
    species_table_1col = []

    # Define the output path
    for specie in species_table:
        if not specie.startswith("#"):
            specie_sp = specie.split("\t")
            species_table_1col.append(specie_sp[0])
    foldername = "SynNetBuild" + datetime.datetime.now().strftime("%Y%m%d_%H%M")
    result_path = f"{output_path}/{foldername}-SynNet-k{hits}s{anchors}m{gaps}"
    os.makedirs(result_path)
    diamond_flodername = "DiamondGenomesDB" + datetime.datetime.now().strftime("%Y%m%d_%H%M")
    diamond_db_path = f"{output_path}/{diamond_flodername}"
    os.makedirs(diamond_db_path)
    os.chdir(result_path)

    # Build a database with diamond
    for sp in species_table_1col:
        print(f"make database for species_{sp}")
        subprocess.run(["diamond", "makedb", "--in", f"{data_path}/{sp}.pep", "-d",
                        f"{diamond_db_path}/{sp}", "-p", str(threads)], check=True)

    # Species-to-species pairwise diamond alignment
    for i in species_table_1col:
        for j in species_table_1col:
            print(f"blast {i} against {j}")
            subprocess.run(
                ["diamond", "blastp", "-q", f"{data_path}/{i}.pep", "-d",
                 f"{diamond_db_path}/{j}", "-o", f"{i}_{j}", "-p",
                 str(threads), "--max-hsps", "1", "-k", str(hits)], check=True)

    # Calculate the synteny blocks between species
    print("Now we do intraspecies MCScanX, prepare inputs.")
    for i in species_table_1col:
        for j in species_table_1col:
            if i == j:
                print(f"{i}_{i} is {i}.blast")
                os.rename(f"{i}_{i}", f"{i}.blast")
                print(f"{i}.bed is {i}.gff")
                shutil.copy(f"{data_path}/{i}.bed", f"{result_path}/{i}.gff")
                print(f"Intra-species MCScanX here for species {i}_{j}")
                subprocess.run(["MCScanX", i, "-s", str(anchors), "-m", str(gaps)], check=True)
                if duplicate is True:
                    subprocess.run(["duplicate_gene_classifier", i], check=True)
                    print("run duplicate_gene_classifier")

                if tandem is True:
                    subprocess.run(
                        ["detect_collinear_tandem_arrays", "-g", f"{result_path}/{i}.gff", "-b",
                         f"{result_path}/{i}.blast", "-c", f"{result_path}/{i}.collinearity", "-o",
                         f"{result_path}/{i}.tandem.collinear"], check=False)
                    print("detect collinear tandem arrays.")

    # Create a non-redundant set
    combinations = list(itertools.combinations(species_table_1col, 2))
    for combo in combinations:
        i = combo[0]
        j = combo[1]
        print(f"merge {i}_{j} {j}_{j} into {i}_{j}.blast")
        merge_files(f"{i}_{j}.blast", [f"{i}_{j}", f"{j}_{i}"])
        print(f"merge {i}.bed {j}.bed into {i}_{j}.gff")
        merge_files(f"{i}_{j}.gff", [f"{data_path}/{i}.bed", f"{data_path}/{j}.bed"])
        print(f"Intra-species MCScanX here for species {i}_{j}")
        subprocess.run(["MCScanX", "-a", "-b 2", f"{i}_{j}", "-s", str(anchors), "-m", str(gaps)], check=True)

    # Output the results
    lines_list = []
    for filename in os.listdir('.'):
        if filename.endswith('.collinearity'):
            with open(filename, 'r') as infile:
                for line in infile:
                    if not line.startswith("#"):
                        lines_list.append(f"{filename}:{line.strip()}")
                    else:
                        lines_list.append(line.rstrip("\n"))

    output_filename = f"SynNet-k{hits}s{anchors}m{gaps}"
    with open(f"{result_path}/{output_filename}", 'a+', encoding='utf-8') as f:
        block_score_table = []
        for line in lines_list:
            if line.startswith("## Alignment"):
                block_score_table.clear()
                part1 = line.split("=")
                part2 = part1[1].split(" ")
                block_score_table.append(part2[0])
            else:
                if "#" not in line:
                    line_sp = line.split("\t")
                    block_id = line_sp[0].replace(".collinearity:", "").replace(" ", "").replace(":", "b")
                    block_score = block_score_table[-1]
                    gene1 = line_sp[1]
                    gene2 = line_sp[2]
                    out_str = f"{block_id}\t{block_score}\t{gene1}\t{gene2}"
                    print(out_str, file=f)


def main():
    global duplicate, tandem
    text1 = '''This Script is for Constructing Synteny Network Database; 
       Usage: python synetfind.py -i species_namelist_file -k 6 -s 5 -m 25 -p 8
       '''
    text2 = '''-i: list of species names; 
       input a list of file names for each species' pep and bed files
       '''
    text3 = '''-d: directory of pep_file and bed_file; 
       Directory where the pep file and bed file are located
       '''
    text4 = '''-o: output_dir; 
       Output files to the specified directory
       '''
    text5 = '''-k: parameter for Diamond: # of top hits (default=6); 
       -k 0: All hits; -k 25: Diamond Default; -k 6: MCScanX suggested
       '''
    text6 = '''-s: parameter for MCScanX: Minimum # of Anchors for a synteny block (default=5); 
       Higher, stricter 
       '''
    text7 = '''-m: parameter for MCScanX: Maximum # of Genes allowed as the GAP between Anchors (default=25); 
       Fewer, stricter
       '''
    text8 = '''-p: Number of threads: Used for Diamond makedb|Diamond blastp|MCSCanX (default=8))
       '''
    parser = argparse.ArgumentParser(description=text1)
    parser.add_argument('-i', '--species_namelist_file', type=str, required=True, help=text2)
    parser.add_argument('-d', '--data_path', type=str, required=True, help=text3)
    parser.add_argument('-o', '--out_dir', type=str, required=True, help=text4)
    parser.add_argument('-k', '--hits', type=int, default=6, help=text5)
    parser.add_argument('-s', '--anchors', type=int, default=5, help=text6)
    parser.add_argument('-m', '--gaps', type=int, default=25, help=text7)
    parser.add_argument('-p', '--threads', type=int, default=8, help=text8)
    parser.add_argument('-D', '--duplicate', action='store_true', help='-D: MCScanX duplicate_gene_classifier')
    parser.add_argument('-T', '--tandem', action='store_true', help='-T: MCScanX detect_collinear_tandem_arrays')
    args = parser.parse_args()
    para = " ".join(sys.argv)
    feedback = f'Execution parameters:\t{para}'
    print(feedback, end="\n", flush=True)
    if args.duplicate:
        duplicate = True
    else:
        duplicate = False
    if args.tandem:
        tandem = True
    else:
        tandem = False

    synet_build(complete_path(args.species_namelist_file), complete_path(args.data_path), complete_path(args.out_dir),
                args.hits, args.anchors, args.gaps, args.threads, duplicate, tandem)


if __name__ == "__main__":
    main()
