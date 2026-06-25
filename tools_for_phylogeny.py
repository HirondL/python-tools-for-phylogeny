import re
import csv
import hashlib
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

def load_from_txt(filepath):
    """
    Load fasta data from a .txt file and return it as a dictionary {"seq_name", "seq", ...}.
    """
    with filepath.open('r') as f:
        L = [line.rstrip() for line in f]

    dict_data = {}
    for line in L:
        if line[0] == ">":
            current_seq = line

            if current_seq in dict_data: # Safety to dodge having duplicate of sequence names
                print("Warning: sequence already loaded", current_seq)

            dict_data[current_seq] = ""
        else:
            if line[-1] == "*":
                dict_data[current_seq] += line[:-1]
            else:
                dict_data[current_seq] += line
    
    return dict_data

def write_to_txt(data, filepath, chunks=True):
    """
    Write in a .txt file fasta data given as a dictionary {"seq_name", "seq", ...}.
    chunks option set to True will write the sequences on multiple lines, otherwise on a single line.
    """
    with filepath.open('w') as f:
        for seq_name in data:
            f.write(seq_name+'\n')
            if not chunks:
                f.write(data[seq_name]+'\n')
            else:
                prot_seq = data[seq_name]
                if len(prot_seq) > 80:
                    for i in range(len(prot_seq)//80): # I took the chunk length used by NCBI
                        f.write(prot_seq[i*80:(i+1)*80]+'\n')
                    last_chunk = prot_seq[(i+1)*80:]
                    if last_chunk != '':
                        f.write(last_chunk+'\n')
                else:
                    f.write(prot_seq+'\n')

def find_index_of_motif(seq, motif):
    """
    For a given sequence seq (str), looks for a motif (str) in it. Tolerant to undetermine 'X' residue and cut '-'.
    If found, return the index (int) where the first matching sequence was found, else False.
    """
    index_seq = 0
    index_motif = 0
    lenght_seq = len(seq)
    lenght_motif = len(motif)
    nb_cuts = 0
    L_index_list = []
    while index_seq < lenght_seq:
        if (seq[index_seq] == motif[index_motif]) or (motif[index_motif] == 'X'):
            index_motif += 1
            if index_motif == lenght_motif:
                index_hit = index_seq - lenght_motif - nb_cuts +1 # to match the indexes
                print(seq[index_hit - 3: index_hit + lenght_motif + 3]) # small window to see the context of the hit
                L_index_list.append(index_hit)
                index_motif = 0
        elif seq[index_seq] == '-':
            nb_cuts += 1
        else:
            nb_cuts = 0
            index_motif = 0
        index_seq += 1
    return L_index_list if L_index_list!=[] else False

def check_motif_at_index(seq, index, motif, cut=False):
    """
    For a given sequence seq (str), check if there is the motif (str) at the index (int).
    Change cut to True to allow cuts tolerance
    Tolerant to undetermine 'X' residue and cut '-' if specified so.
    Return True if the motif is found, else False.
    """
    index_seq = index
    index_motif = 0
    lenght_seq = len(seq)
    lenght_motif = len(motif)
    while index_seq < lenght_seq:
        if (seq[index_seq] == motif[index_motif]) or ((motif[index_motif] == 'X') and (seq[index_seq] != '-')):
            index_motif += 1
            index_seq += 1
            if index_motif == lenght_motif:
                return True
        elif cut and (seq[index_seq] == '-'):
            index_seq += 1
        else:
            return False
    return False #fallback

def write_selected_sequence_to_txt(source_path, output_path, list_selection):
    """
    source_path: path to a .txt file with the data as fasta
    output_path: path to the .txt file the selected sequence will be written in
    list_selection: list of strings with the name of the sequences to select
    """
    with source_path.open('r') as f:
        L = [line.rstrip() for line in f]

    with output_path.open('w') as f:
        in_selected_seq = False
        for line in L:
            if line[0] == ">":
                if line in list_selection:
                    in_selected_seq = True
                    f.write(line+'\n')
                else:
                    in_selected_seq = False
            elif in_selected_seq:
                f.write(line+'\n')

def reformat_NCBI_name(ncbi_name):
    """
    Small tool to parse NCBI fasta sequence description and return groupe and species from the description.
    """
    name = ncbi_name.split(" ")
    prot, gp, sp = name[0][1:], '', ''
    for word in range(len(name)):
        if name[word][0] == '[':
            # A lot of magic numbers based on the format of the name's hit
            gp = name[word][1:]
            sp = name[word+1][:-1] if (name[word+1][-1] == "]" or name[word+1][-1] == ".") else name[word+1]
            break
    return gp, sp, prot

def rename_NCBI_sequence(source_path, output_path, details=''):
    """
    Rename the sequence in source_path file to have them in a format 
        "Group_speciesdetails_protnb NCBI" more suitable to build trees.
    Write a new sequence document in output_path as a .txt
    """
    with source_path.open('r') as f:
        L = [line.rstrip() for line in f]

    with output_path.open('w') as f:
        for line in L:
            if line[0] == ">":
                gp, sp, prot = reformat_NCBI_name(line)
                new_name = '>' + gp + '_' + sp + details + '_' + prot + ' NCBI'
                f.write(new_name+'\n')
            else:
                f.write(line+'\n')

def rename_Phycocosm_sequence(hits_txt, hits_csv, output_path, details=''):
    """
    Rename the sequence in source_path file to have them in a format 
        "Group_speciesdetails_protnb Phycocosm" more suitable to build trees.
    If the name is not found in the hits table, it will print it and keep it the same.
    Write a new sequence document in output_path as a .txt
    """
    with hits_txt.open('r+') as f:
        L = [line.rstrip() for line in f]

    dict_name = {}
    with hits_csv.open('r') as f:
        csvFile = csv.reader(f)
        for lines in csvFile:
            if lines[0] != "Hit":
                # A lot of magic numbers coming for the hits table: 14 --> organism names; 11 --> hits names; 0 --> hits refs
                L_name = lines[14].split(" ")
                gp = L_name[0]
                sp  = L_name[1][:-1] if (L_name[1][-1] == ".") else L_name[1]
                dict_name['>'+lines[11]] = ['>' + gp + '_' + sp, '_' + lines[0] + ' Phycocosm'] # as a list so details could be added in between

    with output_path.open('w') as f:
        for line in L:
            if line[0] == ">":
                if line not in dict_name:
                    print(line)
                    f.write(line+'\n')
                else:
                    new_name = dict_name[line][0] + details + dict_name[line][1]
                    f.write(new_name+'\n')
            else:
                if line[-1] == "*":
                    line = line[:-1]
                if line != '':
                    f.write(line+'\n')


def reformat_Uniprot_name(ncbi_name):
    """
    Small tool to parse NCBI fasta sequence description and return groupe and species from the description.
    """
    name = ncbi_name.split(" ")
    prot, gp, sp = name[0].split('|')[1], '', ''
    for word in range(len(name)):
        if name[word][0:2] == 'OS':
            # A lot of magic numbers based on the format of the name's hit
            gp = name[word][3:]
            sp = name[word+1][:-1] if (name[word+1][-1] == "]" or name[word+1][-1] == ".") else name[word+1]
            break
    return gp, sp, prot

def rename_Uniprot_sequence(source_path, output_path, details=''):
    """
    Rename the sequence in source_path file to have them in a format 
        "Group_speciesdetails_protnb Uniprot" more suitable to build trees.
    Write a new sequence document in output_path as a .txt
    """
    with source_path.open('r') as f:
        L = [line.rstrip() for line in f]

    with output_path.open('w') as f:
        for line in L:
            if line[0] == ">":
                gp, sp, prot = reformat_Uniprot_name(line)
                new_name = '>' + gp + '_' + sp + details + '_' + prot + ' Uniprot'
                f.write(new_name+'\n')
            else:
                f.write(line+'\n')

def reformat_Phaeoexplorer_sequence(source_path, output_path, details=''):
    """
    Reformat the sequence in source_path file to have them in a format 
        "G_speciesdetails_contig Phaeoexplorer" more suitable to build trees.
    Write a new sequence document in output_path as a .txt
    """
    with source_path.open('r') as f:
        L = [line.rstrip() for line in f]

    with output_path.open('w') as f:
        for line in L:
            if line[0] == ">":
                temp1 = line.split('_')
                temp2 = temp1[1].split('-')
                gp, sp = temp2[0], temp2[1]
                contig = temp1[-1].split(' ')[0]
                new_name = '>' + gp + '_' + sp + details + '_' + contig + ' Phaeoexplorer'
                f.write(new_name+'\n')
            else:
                f.write(line+'\n')

def add_details(name, newdetails):
    """
    Given a sequence name in the format "Group_speciesdetails_protnb NCBI", it will add new details (newdetails, str) before the protnb such as:
    "Group_speciesdetails_newdetails_protnb NCBI", return this new name as a string
    """
    new_name = name.split('_')
    new_name[1] = new_name[1] + '_' + newdetails
    return '_'.join(new_name)

def add_details_to_file(in_path, out_path, newdetails):
    """
    Given a file with sequence names in the format "Group_speciesdetails_protnb NCBI", it will add new details (newdetails, str) before the protnb such as:
    "Group_speciesdetails_newdetails_protnb NCBI", return this new name as a string
    Write the renamed sequences in output_path as a .txt
    """
    with in_path.open('r') as f:
        L = [line.rstrip() for line in f]

    with out_path.open('w') as f:
        for line in L:
            if line[0] == ">":
                new_name =add_details(line, newdetails)
                f.write(new_name+'\n')
            else:
                f.write(line+'\n')

def sort_seq_by_tree(dir_path, seq_file, tree_file, chunks=True):
    """
    Given the sequences in 'seq_file' that were used to build a tree in 'tree_file',
    it will write a new sequence file with sequences' order matching the one found in the tree.
    Sequence name format is expected to be: "Group_species_details__protnb database ..."
    chunks option set to True will write the sequences on multiple lines, otherwise on a single line.
    """
    tree_path = dir_path / tree_file
    with tree_path.open('r') as f:
        tree = [line.rstrip() for line in f][0] # when training on SeaView trees, it has only one line
    tree = re.split(r'[:(,]', tree)
    tree = [word for word in tree if len(word)] # to remove empty words
    L_name = [word for word in tree if (word[0].isalpha() or word[0]=='*' or word[0]=='°')] # check if first charcter is a letter or a '*' --> if so, it's a sequence name

    seq_path = dir_path / seq_file
    data = load_from_txt(seq_path)
    dict_seq_to_tree = {seq_name.split()[0][1:]: seq_name for seq_name in data}
    dict_tree_to_seq = {seq_name: seq_name.split()[0][1:] for seq_name in data}

    if (len(data) != len(L_name)): # Safety, extended check
        print('WARNING: Number of items in sequence and tree are not matching.')
        print('Sequence:', len(data))
        print('Tree:', len(L_name))
        print()
        print('In seq but not in tree:')
        for seq in data:
            if dict_tree_to_seq[seq] not in L_name:
                print(dict_tree_to_seq[seq])
        print()
        print('In tree but not in seq:')
        for seq in L_name:
            if seq not in dict_seq_to_tree:
                print(seq)
        print()
        print('Process killed')
        return

    # sligth modification of write_to_txt()
    file_name = seq_file.split('.')[0] + '_sorted.txt'
    file_path = dir_path / file_name
    with file_path.open('w') as f:
        for name in L_name:
            seq_name = dict_seq_to_tree[name]
            f.write(seq_name+'\n')
            if not chunks:
                f.write(data[seq_name]+'\n')
            else:
                prot_seq = data[seq_name]
                if len(prot_seq) > 80:
                    for i in range(len(prot_seq)//80): # I took the chunk length used by NCBI
                        f.write(prot_seq[i*80:(i+1)*80]+'\n')
                    last_chunk = prot_seq[(i+1)*80:]
                    if last_chunk != '':
                        f.write(last_chunk+'\n')
                else:
                    f.write(prot_seq+'\n')

def hash_seq(seq):
    """
    Return the hash for a given str sequence
    """
    return hashlib.md5(seq.encode()).hexdigest()

def compare_datasets_exact(source_data, sink_data, verbose=False):
    """
    Compare an input dataset of sequences that one may want to add to an output dataset of sequences.
    It returns the identifiers of the sequences of the input dataset that are not in the output dataset.
    Source and sink because the idea is to take the ones that are only in the source dataset to add them to the sink dataset.
    Inputs:
        source_data = {seqID: seq(str)}
        sink_data = {seqID: seq(str)}
    Output:
        return the list of sequences' identifiers not matching the other dataset as the first element
        then the list of the identifiers of sequences that were matching the other dataset as the second element
    """
    source_hash = {seq: hash_seq(source_data[seq]) for seq in source_data}
    sink_hash = {seq: hash_seq(sink_data[seq]) for seq in sink_data}
    new_seq = [seq for seq in source_hash if source_hash[seq] not in list(sink_hash.values())]
    known_seq = [seq for seq in source_hash if seq not in new_seq]
    if verbose:
        # For each known sequence, it will print the name of the sequence matching in sink_data.
        # Thus, if the name in the source_data is better, one can see it and decide to keep it instead of the other one.
        for seq in known_seq:
            print(seq, [sink_seq for sink_seq in sink_hash if sink_hash[sink_seq]==source_hash[seq]])
    return new_seq, known_seq

def compare_dataset_identity(data, cuts=False, threshold=0.99, get_matrix=False, verbose=False, visual=False):
    """
    Compare the sequences in a dataset and compute a pairwise identity score.
    It returns a list of scores with the associated sequences' names above a given threshold.
    Inputs:
        data = {seqID: seq(str)}, sequences aligned
        cuts: set to True if one wants to count 2 '-' as identical between 2 sequences
        threshold (float): pair of sequences with an identity score above this threshold will be consider as postive
        verbose: set to True if one wants the positive pair of sequences to be dipslay with their identity scores
        visual: set to True if one wants to display the identity matrix.
            Note that the score matrix is symmetrical and have 1s on the diagonal. Here I chose to display only the upper part.
            If you want a proper score matrix, get the score matrix M from this function and do M + Id + M.T
        get_matrix: set to True to export the score matrix as well
    Output:
        if get_matrix==False:
            return the positive list: [[score, seq_id_1, seq_id_2], [...], ...]
        if get_matrix==True:
            return the postive list plus the score matrix: [[score, seq_id_1, seq_id_2], [...], ...], score_matrix (numpy array)
    """
    #---------------Building the score matrix---------------
    size = len(data)
    score_matrix = np.zeros((size,size))
    seq_id = list(data.keys())
    len_seq = len(data[seq_id[0]])

    for seq_base in range(size):
        for seq_test in range(seq_base+1, size):
            if cuts:
                score_matrix[seq_base,seq_test] = sum(x == y for x, y in zip(data[seq_id[seq_base]], data[seq_id[seq_test]]))/len_seq
            else:
                similarity_nb = sum((x == y and (x!='-' and y!='-')) for x, y in zip(data[seq_id[seq_base]], data[seq_id[seq_test]]))
                non_double_cuts_nb = sum(not(x=='-' and y=='-') for x, y in zip(data[seq_id[seq_base]], data[seq_id[seq_test]]))
                score_matrix[seq_base,seq_test] = similarity_nb / non_double_cuts_nb
    
    #---------------Displaying the score matrix---------------
    if visual:
        fig, ax = plt.subplots()
        mat = ax.matshow(score_matrix, cmap='Spectral', vmin=0, vmax=1)
        ax.set_title("Sequences normalized identity scores")
        fig.colorbar(mat, orientation='vertical')
        plt.tight_layout()
    
    #---------------Keeping the sequences with a score above the threshold and printing their IDs---------------
    threshold_matrix = ((score_matrix >= threshold) * score_matrix) # Numpy trick to keep only value >= 0.99
    coordinates = np.nonzero(threshold_matrix)

    if verbose:
        print('Threshold >=', threshold)
    
    L_high_similarity = []
    for i in range(len(coordinates[0])):
        current_match = [threshold_matrix[coordinates[0][i],coordinates[1][i]], seq_id[coordinates[0][i]],seq_id[coordinates[1][i]]]
        L_high_similarity.append(current_match)
        if verbose:
            print(current_match[0], current_match[1], current_match[2])
    
    #---------------OUTPUT---------------
    if get_matrix:
        return L_high_similarity, score_matrix
    return L_high_similarity
    

def spot_mismatch(seq_1, seq_2, verbose=False):
    """
    For 2 given aligned sequences, return a dictionary of the found mismatch: {index: [aa_1, aa_2], ...}
    Set verbose to True to directly print the dictionary.
    Will print an alert if more than 20 mismatch are found.
    """
    dict_mismatch = {}
    len_seq = len(seq_1)

    for aa in range(len_seq):
        if seq_1[aa] != seq_2[aa]:
            dict_mismatch[str(aa)] = [seq_1[aa], seq_2[aa]]
    
    if len(dict_mismatch) > 20:
        print('WARNING: more than 20 mismatches found!')
    
    if verbose:
        print('Found '+str(len(dict_mismatch))+' mismatch(es).')
        for mismatch in dict_mismatch:
            print(mismatch+':', dict_mismatch[mismatch][0], '<->', dict_mismatch[mismatch][1])
    
    return dict_mismatch


def plot_sequence_stats(dict_list, labels=None):
    """
    Plot ECDF and density estimation for sequence lengths in each dictionary.

    Parameters:
    - dict_list: List of dictionaries, each containing sequences.
    - labels: Optional list of labels for each dictionary (for legend).
    """
    if labels is None:
        labels = [f"Dict {i+1}" for i in range(len(dict_list))]

    # Collect all sequence lengths for each dictionary
    all_lengths = []
    for d in dict_list:
        lengths = [len(seq) for seq in d.values()]
        all_lengths.append(lengths)

    # Create figure and subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Plot ECDF
    for lengths, label in zip(all_lengths, labels):
        sorted_lengths = np.sort(lengths)
        ecdf = np.arange(1, len(sorted_lengths)+1) / len(sorted_lengths)
        ax1.step(sorted_lengths, ecdf, label=label, where='post')
    ax1.set_title("ECDF of Sequence Lengths")
    ax1.set_xlabel("Sequence Length")
    ax1.set_ylabel("ECDF")
    ax1.legend()
    ax1.grid(True)

    # Plot Density Estimation
    for lengths, label in zip(all_lengths, labels):
        kde = gaussian_kde(lengths)
        x_grid = np.linspace(min(lengths)-1, max(lengths)+1, 1000)
        ax2.plot(x_grid, kde(x_grid), label=label)
    ax2.set_title("Density Estimation of Sequence Lengths")
    ax2.set_xlabel("Sequence Length")
    ax2.set_ylabel("Density")
    ax2.legend()
    ax2.grid(True)

    plt.tight_layout()
    plt.show()