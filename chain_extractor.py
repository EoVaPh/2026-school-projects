import argparse

parser = argparse.ArgumentParser(
    prog='Chain extractor',
    description='''Extract chains or their fragments from PDB files.''',
    epilog='(C) 2026 Egor Vasilenko')

parser.add_argument('-f', '--file', type=str)
parser.add_argument('-c', '--chain', type=str)
parser.add_argument('-b', '--begin', type=int, default=None)
parser.add_argument('-e', '--end', type=int, default=None)

args = parser.parse_args()

file = args.file
chain = args.chain
begin = args.begin
end = args.end

structure_file = open(file, 'r')
structure_lines = structure_file.readlines()
structure_file.close()

for line in structure_lines:
    if 'ATOM' in line and ' ' + chain + ' ' in line:
        res_id = int(line.split()[5])
        print_flag = True
        if begin is not None and res_id < begin:
            print_flag = False
        if end is not None and res_id > end:
            print_flag = False
        if print_flag:
            print(line, end='')
