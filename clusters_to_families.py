def read_structure_protein_name(structure_file_name: str) -> str:
    '''Find line in a structure file that probably contains the name of a
       polypeptide.'''

    structure_file = open('amyloid_structures/' + structure_file_name, 'r')
    structure_lines = structure_file.readlines()
    structure_file.close()

    for line in structure_lines:
        if '.pdb' in structure_file_name and 'TITLE' in line:
            return line.replace('TITLE', '').strip()

        if '.cif' in structure_file_name and \
           ('_entity_name_com.name' in line or \
            '_entity.pdbx_description' in line):
            return line.replace('_entity_name_com.name', '').\
                        replace('_entity.pdbx_description', '').\
                        strip()


clusters_file = open('clustered_sequences.txt.clstr', 'r')
clusters_lines = clusters_file.readlines()
clusters_file.close()

families_file = open('families.fa', 'w')

for line in clusters_lines:
    if line.strip()[0] == '>':
        family_header_flag = True
    else:
        structure_file_name = line.strip().split('>')[1][:8]
        if family_header_flag:
            families_file.write(
                '>' + read_structure_protein_name(structure_file_name) + '\n'
            )
        family_header_flag = False
        families_file.write(structure_file_name[:-4] + '\n')

families_file.close()
