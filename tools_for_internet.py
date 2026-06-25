import json
import subprocess
import tools_for_phylogeny as tools
from pathlib import Path # not sure if I need it to deal with path given as inputs here...

from Bio import Entrez, SeqIO
from urllib.error import HTTPError
Entrez.email = "your@email.com"  # required by NCBI

def UniProt_call_database(prot_id):
    web_call = subprocess.run(['curl', 'https://rest.uniprot.org/uniprotkb/'+prot_id], capture_output=True)
    # Convert bytes → string
    response_text = web_call.stdout.decode('utf-8')
    # Parse JSON → Python dict
    return json.loads(response_text)

def UniProt_get_prot_info(prot_id):
    data = UniProt_call_database(prot_id)
    if 'sequence' in data:
        name_seq = '>' + '_'.join(data['organism']['scientificName'].split(' ')[:2]) + '_' + data['primaryAccession'] + ' UniProt'
        seq = data['sequence']['value']
        return [True, name_seq, seq]
    return [False, prot_id, data]

def UniProt_get_all_prot(L_prot_ids):
    dict_prot = {}
    dict_failure = {}
    for prot_id in L_prot_ids:
        result = UniProt_get_prot_info(prot_id)
        if result[0]:
            dict_prot[result[1]] = result[2]
        else:
            dict_failure[prot_id] = result[2]      
    return dict_prot, dict_failure

def NCBI_call_database(prot_id, base='protein'):
    try:
        handle = Entrez.efetch(
            db=base, # can also use 'nucleotide' to get DNA sequences.
            id=prot_id,
            rettype="fasta",
            retmode="text"
        )
        record = SeqIO.read(handle, "fasta-blast")
    except (HTTPError, ValueError):
        record = False  # no result

    return record

def NCBI_get_prot_info(prot_id, base='protein'): # can also use 'nucleotide' to get DNA sequences.
    data = NCBI_call_database(prot_id, base=base)
    if data:
        if base=='protein':
            gp, sp, prot = tools.reformat_NCBI_name('>'+data.description)
            name_seq = '>' + gp + '_' + sp + '_' + prot + ' NCBI'
        else:
            name_seq = '>'+data.description
        seq = str(data.seq)
        return [True, name_seq, seq]
    return [False, prot_id, 'Not found in NCBI using the Entrez.efectch from biopython.']

def NCBI_get_all_prot(L_prot_ids, base='protein'): # can also use 'nucleotide' to get DNA sequences.
    dict_prot = {}
    dict_failure = {}
    for prot_id in L_prot_ids:
        result = NCBI_get_prot_info(prot_id, base=base)
        if result[0] == False:
            dict_failure[prot_id] = result[2]
        else:
            dict_prot[result[1]] = result[2]
    return dict_prot, dict_failure

def ATGC_encode(dirpath, filename, name_threshold=30, name_sep='@', name_tag='!'):
    """
    The rationale here is to get a tree computed with ATGC Montpellier cluster. But it uses Phylip format that trims names bigger that 30 caracters.
    Thus, given a sequence TXT file in FASTA format with the format used in this work (>Group_species_details_protnb):
        1) It will create an identical TXT file except with '_short' in the name of the file, and with sequence names as >@!Gsp_details_protnb@.
        2) It will also create a JSON file containing the dictionary {short_name: long_name} so one can decode the tree after computation:
            Split the tree using (eg) '|', then look for 1st caracter (eg) '!', then replace with the corresponding long name.
    Rk: '#' seems to cause problem on the ATGC cluster. Careful with '|' if using sequences from Phycocosm.
    Rk: no check to verify that all names are different: protnb should ensure that.
    """
    # Load the intial data
    file_path = dirpath / filename
    data = tools.load_from_txt(file_path)

    # Shorten Group_species_details_protnb into Gsp_details_protnb
    L_names = list(data.keys())
    dict_trimmed_names = {names: names.split(' ')[0] for names in L_names}
    dict_names = {dict_trimmed_names[names]: dict_trimmed_names[names].split('_') for names in dict_trimmed_names} # splits to get [Group, species, details, protnb]
    convert_dict_short_names = {} # dict long to short names
    print('Please consider renaming the following sequences:')
    for name in dict_names:
        L_short = dict_names[name]
        if L_short[0][1] == '*':
            short_name = '>'+name_sep+name_tag + L_short[0][1:3] + L_short[1][:2] + '_'.join(L_short[2:]) + name_sep
        else:
            short_name = '>'+name_sep+name_tag + L_short[0][1:2] + L_short[1][:2] + '_'.join(L_short[2:]) + name_sep
        if len(short_name) > (name_threshold + 1): # because you have this special caracter '>' that won't be in the Phylip format
            print(short_name)
        convert_dict_short_names[name] = short_name
    reverse_dict_short_names = {convert_dict_short_names[name]: name for name in convert_dict_short_names} # dict short to long names

    # Export the dictionary used to reverse the process as a JSON file
    reverse_filename = 'reverse_' + filename.split('.')[0] + '_short.json'
    reverse_path = dirpath / reverse_filename
    with open(reverse_path, "w") as f:
            json.dump(reverse_dict_short_names, f)

    # Export the sequences data with short names
    short_data = {convert_dict_short_names[dict_trimmed_names[name]]: data[name] for name in dict_trimmed_names} # sequences data with short names
    short_filename = filename.split('.')[0] + '_short.txt'
    out_path = dirpath / short_filename
    tools.write_to_txt(short_data, out_path)

def ATGC_decode(dirpath, treename, jsonname, name_sep='@', name_tag='!'):
    """
    So one got the tree computed with ATGC Montpellier cluster. Because it uses Phylip format that trims names bigger that 30 caracters,
    one used ATGC_encode() just above to have short sequence names.
    Now, giving back the TXT tree and JSON decoding dictionary, this function will split the tree using (eg) '|',
    then look for 1st caracter (eg) '!', then replace with the corresponding long name to get back >Group_species_details_protnb.
    Finally, it will write the new TXT tree as '..._renamed.txt'.
    """
    # Load the tree
    tree_path = dirpath / treename
    with tree_path.open('r') as f:
        tree = [line.rstrip() for line in f][0] # when training on SeaView trees, it has only one line, maybe it has more on ATGC cluster...
    
    # Load the decoding dictionary
    reverse_path = dirpath / jsonname
    with open(reverse_path, 'r') as file:
        reverse_dict = json.load(file)

    # Replace short names by long names
    split_tree = tree.split(name_sep)
    new_tree = []
    for i in range(len(split_tree)):
        if split_tree[i][0] == name_tag:
            short_name = '>'+name_sep + split_tree[i] + name_sep
            new_tree.append(reverse_dict[short_name][1:])
        else:
            new_tree.append(split_tree[i])
    new_tree = ''.join(new_tree)

    # Export the renamed tree
    new_treename = treename.split('.')[0] + '_renamed.txt'
    out_path = dirpath / new_treename
    with out_path.open('w') as f:
        f.write(new_tree)
