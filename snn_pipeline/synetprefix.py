# CEG
# haochenli
# Python3
# date: 24/7/2024 上午11:05
import argparse
import os
import datetime


def complete_path(input_path):
    if not os.path.isabs(input_path):
        input_path = os.path.join(os.getcwd(), input_path)
    return input_path


def read_table(file_name):
    out_table = []
    in_put = open(file_name, "r")
    for line1 in in_put:
        line1_rs1 = line1.rstrip("\n")
        line1_rs2 = line1_rs1.strip(" ")
        out_table.append(line1_rs2)
    in_put.close()
    return out_table


def block_trans(input_file, output_dir):
    foldername = "SynNetTrans" + datetime.datetime.now().strftime("%Y%m%d_%H%M")
    os.makedirs(f"{output_dir}/{foldername}")
    input_file_name = os.path.basename(input_file)
    print("Creating prefix file", end="\n", flush=True)
    with open(f"{output_dir}/{foldername}/{input_file_name}.prefix", 'a+', encoding='utf-8') as f:
        merged_block_table = read_table(input_file)
        block_id_table = []
        r = 0
        for l1 in merged_block_table:
            r = r + 1
            l1_sp = l1.split("\t")
            col1_sp = l1_sp[0].split("-")
            block_id = col1_sp[0]
            if block_id not in block_id_table and r == 1:
                block_id_table.append(block_id)
                out_str = block_id + "\t" + "\t".join(l1_sp[2:])
                print(out_str, file=f, end=">")
            elif block_id in block_id_table:
                out_str = "\t".join(l1_sp[2:])
                print(out_str, file=f, end=">")
            else:
                block_id_table.clear()
                block_id_table.append(block_id)
                out_str = "\n" + block_id + "\t" + "\t".join(l1_sp[2:])
                print(out_str, file=f, end=">")


def main():
    text1 = '''The script converts the SynNet_file from a two-column format to a 
    format where collinearity genes of same block are on one line.; 
          Usage: python Merged_block_prefix -n SynNet_file
          '''
    text2 = '''-n: input Synnet_file; 
    input Synnet_file which generated through SynBuild-X script
    '''
    text3 = '''-o: output_dir; 
       Output files to the specified directory
    '''
    parser = argparse.ArgumentParser(description=text1)
    parser.add_argument('-n', '--synnet_file', type=str, required=True, help=text2)
    parser.add_argument('-o', '--out_dir', type=str, required=True, help=text3)
    args = parser.parse_args()
    block_trans(complete_path(args.synnet_file), complete_path(args.out_dir))


if __name__ == "__main__":
    main()
