# Import statement
import math
import os
import pandas as pd
import requests
import ssl
import subprocess
import time
import wget

ssl._create_default_https_context = ssl._create_unverified_context

### Input data - Make sure everything in here is correct else the program will probably not work
# Name of the project
project_name = "Delta_boucle"

# URL to your job containing the password if there is one,
# exemple : "https://yanglab.nankai.edu.cn/mTM-align/output/mTMxxxxxx/xxxxxx/"
url = "https://yanglab.nankai.edu.cn/mTM-align/output/mTM012917/svcekm/"
# Are the results files already downloaded and placed in the correct folder? True or False (make sure the caps is there)
already_dl = False

### OPTIONS
csv_generation = False
see_aa = True  # Do you want to see similar aminoacids in the pymol session? True or False (make sure the caps is there)
# Color and representation scheme
color_identical = "green"  # Color of the identical amino acids
color_same = "black"  # Color of the same class amino acids
bg_color = "white"  # Color of the background
representation_identical = "stick"  # Represetion of the identical amino acids
representation_same = "line"  # Representtion of the same class amino acids
representation_ligand = "stick"  # Representation of the ligands
# If an error occur during the superposition because of time. Usually 2 is enough
wait_multiplier = 2
# Choose the method for coloration from these two :  "RMSD" or "HD" (for highest_distances)
# Attention : RMSD method doesn't work with only two structures
method_coloration = "RMSD"
# Choose the cutoff in ångström for the coloration. Usually 10 for "HD" and 2 for "RMSD"
cutoff = 2
# Path to PyMOL
PyMolpath = r"C:\Users\Julien\PyMOL\PyMOLWin.exe"


# Main function
def main():
    start = time.perf_counter()

    global path_to_fasta, path_to_newpdb, pdb_dict, result, pos

    # Initialization of the project : creation of necessary folders

    print(f"Initialization of project named '{project_name}'")
    project_folder = create_folder(os.getcwd(), project_name)
    structure_folder = create_folder(os.getcwd(), "all_structures")
    fastas_folder = create_folder(project_folder, "alignements")
    pymol_session = create_folder(project_folder, "session_pymol")
    pymol_script_folder = create_folder(project_folder, "pymol_script")
    processed_folder = create_folder(project_folder, "processed_pdbs")
    aligned_structures_folder = create_folder(project_folder, "super_pdbs")
    csv_folder = create_folder(project_folder, "table_files")
    log_folder = create_folder(project_folder, "log_files")

    print()

    # Managing results folder, and download files from the result on MSTA
    if already_dl:
        msta_results_folder = project_folder + "\\MSTA_raw_results"
        path_to_fasta = msta_results_folder + f"\\seq.fasta"
        path_to_newpdb = msta_results_folder + f"\\new.pdb"
    else:
        try:
            msta_results_folder = create_folder_destructive(project_folder, "MSTA_raw_results")
            start_fasta = time.perf_counter()
            path_to_fasta = download_wget(url, "seq.fasta", msta_results_folder)
            finish_fasta = time.perf_counter()
            start_pdb = time.perf_counter()
            path_to_newpdb = download_wget(url, "new.pdb", msta_results_folder)
            finish_pdb = time.perf_counter()
        except FileNotFoundError:
            print("Error, The files are not present, please check the variable already_dl")
            exit(-1)

    pdb_list = get_list_from_ali(path_to_fasta)

    structure_left_behind = download_pdbs(pdb_list, structure_folder)
    for i in structure_left_behind:
        if os.path.isfile(structure_folder + f"\\{i}.pdb"):
            pass
        else:
            print(f"The structure folder is lacking {i}.pdb, "
                  f"please add it in the folder {project_name}\\structure_folder")
    if check_if_done(project_name, structure_folder, pdb_list):
        pass
    print()
    print("Processing structures for seleno-methionines.")
    mse_to_met(structure_folder, processed_folder, pdb_list)
    print("Done")
    print()
    print("Superimposition of selected structures to the raw result of MSTA")
    print("Using PyMOL!")

    with open(pymol_script_folder + f"\\super{project_name}.pml", "w") as output:
        print(f"""reinitialize
                load {path_to_newpdb}""", file=output)
        for x in pdb_list:
            print(f"""load {processed_folder}\\{x}.pdb
            super {x} and chain A, new""", file=output)
        for x in pdb_list:
            print(f"""save {aligned_structures_folder}\\{x}.pdb, {x}""", file=output)
        print("""quit""", file=output)
    print(f"Waiting for {len(pdb_list) * 0.5 * wait_multiplier} seconds, time to superimpose all structures.")

    subprocess.call([PyMolpath, pymol_script_folder + f"\\super{project_name}.pml"])

    time.sleep(len(pdb_list) * 0.5 * wait_multiplier)
    print("Done!")

    for i in pdb_list:
        print(f"Creation of {i}.fasta")
        align_separate_sequences(path_to_fasta, i, fastas_folder)
    print()
    hm_fasta = how_many_pdb(fastas_folder)

    try:
        pdb_dict = dictionary_coordinates_pdb(pdb_list, aligned_structures_folder)
    except FileNotFoundError:
        print(
            "Error! Superpositon was not completed in time for the next step. "
            "Increase wait_multiplier variable in input and retry")
        exit(-1)

    print("Dictionary of structures and coordinates generated correctly")

    matrix_pdb = matrix_pdbs(pdb_list)
    print("Matrices of structures generated correctly")
    align_len = read_lenght(pdb_list[0], fastas_folder)
    print()

    with open(pymol_script_folder + f"\\{project_name}.pml", "w") as f:
        print("""reinitialize""", file=f)
        for n in pdb_list:
            m = n.replace(".pdb", "")
            print(f"""load {aligned_structures_folder}\\{n}.pdb
    extract {m}_A, chain A and {m}
    delete {m}""", file=f)
        print("""color grey40, all
            select water, resn HOH
            remove water""", file=f)

    with open(pymol_script_folder + f"\\{project_name}_aminoacids.pml", "w"):
        pass

    print("Alignement reading and processing")
    plus_grande_valeur = []
    for chara in range(align_len):
        number = 0
        aa_class = set()
        aa_type = set()
        position = chara + 1
        print(f"Reading characters {position} in alignement...")
        for structure in pdb_list:
            if check_pos_ifchar(chara, structure, fastas_folder):
                number += 1
            aa_type.update(give_char(chara, structure, fastas_folder))
            aa_class.update(str(check_pos_whichclass(chara, structure, fastas_folder)))

        if number >= hm_fasta:
            dict_dist = {}
            for tuplex in range(len(matrix_pdb)):
                x = matrix_pdb[tuplex]
                une, deux = x.split(",")

                a = une.strip()
                b = deux.strip()

                # structure 1 - a
                seq_a = fetch_sequence_in_fasta(a, fastas_folder)
                dash_a = count_dash(seq_a, 0, position)
                calculus_a = position - dash_a
                x1, y1, z1, seqres1 = xyzseq_pdb(pdb_dict, a, calculus_a)

                # structure 2 - b
                seq_b = fetch_sequence_in_fasta(b, fastas_folder)
                dash_b = count_dash(seq_b, 0, position)
                calculus_b = position - dash_b
                x2, y2, z2, seqres2 = xyzseq_pdb(pdb_dict, b, calculus_b)

                dict_dist.update({f"{a}, {seqres1}, {b}, {seqres2}": f"{position}, {distance(x1, y1, z1, x2, y2, z2)}"})

            plus_grande_valeur.append(max(dict_dist.values()))

            if method_coloration == "HD":
                value = max(dict_dist.values())
                pos, result = value.split(",")
                pos, result = int(pos.strip()), float(result.strip())

            elif method_coloration == "RMSD":
                value_list = []
                for i in dict_dist:
                    useless, value = dict_dist[i].split(",")
                    value_list.append(value.strip())
                result = rmsd_list(value_list)
                value_2 = max(dict_dist.values())
                pos, dismax = value_2.split(",")
                pos, dismax = int(pos.strip()), float(dismax.strip())

            else:
                print("Method not specified. Stopping...")
                time.sleep(2)
                exit(-1)

            # Addition of wanted lines in the PyMOL script depending of the maximum distance
            with open(pymol_script_folder + f"\\{project_name}.pml", "a") as f:
                with open(pymol_script_folder + f"\\{project_name}_aminoacids.pml", "a") as f2:
                    for i in pdb_list:
                        if result <= 0.1 * cutoff:
                            print_in_script(i, pos, "br0", f, pdb_dict, fastas_folder)
                            continue
                        elif result <= 0.2 * cutoff:
                            print_in_script(i, pos, "br1", f, pdb_dict, fastas_folder)
                            continue
                        elif result <= 0.3 * cutoff:
                            print_in_script(i, pos, "br2", f, pdb_dict, fastas_folder)
                            continue
                        elif result <= 0.4 * cutoff:
                            print_in_script(i, pos, "br3", f, pdb_dict, fastas_folder)
                            continue
                        elif result <= 0.5 * cutoff:
                            print_in_script(i, pos, "br4", f, pdb_dict, fastas_folder)
                            continue
                        elif result <= 0.6 * cutoff:
                            print_in_script(i, pos, "br5", f, pdb_dict, fastas_folder)
                            continue
                        elif result <= 0.7 * cutoff:
                            print_in_script(i, pos, "br6", f, pdb_dict, fastas_folder)
                            continue
                        elif result <= 0.8 * cutoff:
                            print_in_script(i, pos, "br7", f, pdb_dict, fastas_folder)
                            continue
                        elif result <= 0.9 * cutoff:
                            print_in_script(i, pos, "br8", f, pdb_dict, fastas_folder)
                            continue
                        elif result > 0.9 * cutoff:
                            print_in_script(i, pos, "br9", f, pdb_dict, fastas_folder)
                            continue

                    for i in pdb_list:
                        if len(aa_class) == 1 and len(aa_type) == 1:
                            print_in_script_aminoacids(i, position, color_identical, f2, pdb_dict, fastas_folder,
                                                       representation_identical)
                        if len(aa_class) == 1 and len(aa_type) > 1:
                            print_in_script_aminoacids(i, position, color_same, f2, pdb_dict, fastas_folder,
                                                       representation_same)

        else:
            pass

    if csv_generation:
        csv_file_creation(fastas_folder, project_name, csv_folder)
        read_file = pd.read_csv(csv_folder + f"\\{project_name}.csv", sep=",", header=None)
        read_file.to_excel(csv_folder + f"\\{project_name}.xlsx")

    with open(pymol_script_folder + f"\\{project_name}.pml", "a") as f:
        if see_aa:
            print(f"""run {pymol_script_folder}\\{project_name}_aminoacids.pml""", file=f)
        print(f"""set cartoon_flat_sheets, 1
set cartoon_side_chain_helper, 1
set cartoon_fancy_helices, 1
select Identical_Residues, color {color_identical}
select SameClass_Residues, color {color_same}
hide (hetatm)
show {representation_ligand}, hetatm
orient all
bg_color {bg_color}
util.cnc
deselect
save {pymol_session}\\{project_name}.pse""", file=f)

    print()
    print("Opening PyMOL session")

    subprocess.call([PyMolpath, pymol_script_folder + f"\\{project_name}.pml"])

    finish = time.perf_counter()

    with open(log_folder + f"\\options_{project_name}.txt", "w") as logfile:
        print(f"""
Project name = {project_name}
URL to the mTM-align job = {url}
Structure list = {pdb_list}

Identical amino acid coloration = {color_identical}
Identical amino acid representation = {representation_identical}
Same class amino acid color {representation_same}
Same class amino acid representation = {representation_same}
Ligand representation = {representation_ligand}

Coloration method = {method_coloration}
Cutoff aplied = {cutoff}

Script execution time finished in {round(finish - start, 2)} second(s)
Details:
Fasta file download finished in {round(finish_fasta - start_fasta, 2)} second(s)
PDB file download finished in {round(finish_pdb - start_pdb, 2)} second(s)
""", file=logfile)

    print(f"finished in {round(finish - start, 2)} second(s)")


# Auxiliary functions
def rmsd_list(list_val):
    size = len(list_val)
    mean_intermediate = 0
    for i in list_val:
        mean_intermediate += float(i)
    mean = mean_intermediate / size

    calculus = 0
    for i in list_val:
        calculus += ((float(i) - mean) ** 2)

    the_result = math.sqrt(calculus / (size - 1))
    return the_result


def dummy(something):
    return something


def special_for_this(pdb, position, pdb_dictionnary, fasta_folder):
    seq = fetch_sequence_in_fasta(pdb, fasta_folder)
    dash = count_dash(seq, 0, position)
    true_pos = position - dash
    x, y, z, seqres = xyzseq_pdb(pdb_dictionnary, pdb, true_pos)
    return seqres


def print_in_script_aminoacids(pdb, position, color, file, pdb_dictionnary, fasta_folder, representation):
    seq = fetch_sequence_in_fasta(pdb, fasta_folder)
    dash = count_dash(seq, 0, position)
    true_pos = position - dash
    x, y, z, seqres = xyzseq_pdb(pdb_dictionnary, pdb, true_pos)
    i = pdb.replace(".pdb", "")
    print(f"""show {representation}, {i}_A and resi {seqres}
color {color}, {i}_A and resi {seqres} and not (name N+O+C+CA)""", file=file)


def print_in_script(pdb, position, color, file, pdb_dictionnary, fasta_folder):
    seq = fetch_sequence_in_fasta(pdb, fasta_folder)
    dash = count_dash(seq, 0, position)
    true_pos = position - dash
    x, y, z, seqres = xyzseq_pdb(pdb_dictionnary, pdb, true_pos)
    i = pdb.replace(".pdb", "")
    print(f"""color {color}, {i}_A and resi {seqres} and (name c+o+n+ca)""", file=file)


def csv_file_creation(alignement_folder, project, output_folder):
    fasta_files = [name for name in os.listdir(alignement_folder) if
                   os.path.isfile(os.path.join(alignement_folder, name))]
    print()
    print("The CSV file is generating")
    with open(f"{output_folder}\\{project}.csv", "w") as o_f:
        for file in fasta_files:
            with open(f"{alignement_folder}\\{file}") as in_f:
                for line in in_f:
                    if line[:1].strip() == ">":
                        continue
                    else:
                        name_len = len(file.replace(".fasta", ""))
                        chars = [x for x in line.strip().replace(f"{file}.fasta", "")]
                        chars_to_print = ""
                        for i in chars:
                            chars_to_print += f",{i}"
                        print(file[:name_len] + chars_to_print, file=o_f)

    liste_seq = []
    with open(f"{output_folder}\\{project}.csv", "r") as w_f:
        for line in w_f:
            for file in fasta_files:
                if line.startswith(file.replace(".fasta", "")):
                    name_len = len(file.replace(".fasta", ""))
                    liste = [x for x in line[name_len:].strip().replace(",", "")]
                    liste_seq.append(liste)
    print("Starting comparison of amino acids")
    align_len = read_lenght(fasta_files[0].replace(".fasta", ""), alignement_folder)
    with open(f"{output_folder}\\{project_name}.csv", "a") as w_f:
        string = ","
        for chara in range(align_len):
            if chara != align_len - 1:
                aa_set = set()
                class_set = set()
                for liste in liste_seq:
                    aa_set.update(str(liste[chara]))
                    class_set.update(str(which_class(liste[chara])))
                if '-' in aa_set:
                    string += "0,"
                elif '-' not in aa_set and len(class_set) == 1 and len(aa_set) == 1:
                    string += "3,"
                elif '-' not in aa_set and len(class_set) == 1 and len(aa_set) > 1:
                    string += "2,"
                elif '-' not in aa_set:
                    string += "1,"
            else:
                aa_set = set()
                class_set = set()
                for liste in liste_seq:
                    aa_set.update(str(liste[chara]))
                    class_set.update(str(which_class(liste[chara])))
                if '-' in aa_set:
                    string += "0"
                elif '-' not in aa_set and len(class_set) == 1 and len(aa_set) == 1:
                    string += "3"
                elif '-' not in aa_set and len(class_set) == 1 and len(aa_set) > 1:
                    string += "2"
                elif '-' not in aa_set:
                    string += "1"

        print(string.strip(), file=w_f)
    print("CSV file generated successfully")
    return output_folder + f"\\{project_name}.csv"


def which_class(nature):
    aa_1 = ["A", "G"]  # nonpolar_small
    aa_2 = ["V", "L", "I", "M", "P"]  # nonpolar_medium
    aa_3 = ["F", "W", "Y"]  # aromatic
    aa_4 = ["S", "T", "C", "N", "Q"]  # polar
    aa_5 = ["K", "R", "H"]  # plus
    aa_6 = ["D", "E"]  # minus

    if nature in aa_1:
        return 1
    if nature in aa_2:
        return 2
    if nature in aa_3:
        return 3
    if nature in aa_4:
        return 4
    if nature in aa_5:
        return 5
    if nature in aa_6:
        return 6
    else:
        return 0


def mse_to_met(original_folder, processed_folder, liste):
    """change the selenomethionine in methionine for the pdb files in a folder,
    write new files into a destination folder"""
    for file in os.listdir(original_folder):
        if file[:-4] in liste:
            with open(f"{original_folder}\\{file}", "r") as file_input:
                print(f"Processing {file}")
                with open(f"{processed_folder}\\{file}", "w") as file_output:
                    for line in file_input:
                        if line[:6].strip() == "HETATM" and line[17:20].strip() == "MSE":
                            print(line.strip().replace("HETATM", "ATOM  ").replace("MSE", "MET").replace("SE", " S"),
                                  file=file_output)
                        else:
                            print(line.strip(), file=file_output)


def xyzseq_pdb(dic, pdb, number):
    x, y, z, nseq = (dic.get(f"{pdb}, {number}")).split(',')
    return float(x), float(y), float(z), int(nseq)


def distance(x1, y1, z1, x2, y2, z2):
    return math.sqrt(((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2))


def count_dash(string, start, end):
    return string.count("-", int(start), int(end))


def fetch_sequence_in_fasta(file, fasta_folder):
    with open(f"{fasta_folder}\\{file}.fasta") as f:
        for line in f:
            if line[:1].strip() != ">":
                return line.strip()


def is_aa(sequence, position):
    aminoacids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                  'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    if sequence[position].strip() in aminoacids:
        return True
    else:
        return False


def check_pos_ifchar(number, file, fasta_folder):
    aminoacids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                  'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    with open(f"{fasta_folder}\\{file}.fasta") as f:
        for line in f:
            if line[:1].strip() != ">":
                if line[number].strip() in aminoacids:
                    return True
                else:
                    return False


def give_char(number, file, fasta_folder):
    with open(f"{fasta_folder}\\{file}.fasta") as f:
        for line in f:
            if line[:1].strip() != ">":
                return line[number].strip()


def check_pos_whichclass(number, file, fasta_folder):
    aa_1 = ["A", "G"]  # nonpolar_small
    aa_2 = ["V", "L", "I", "M", "P"]  # nonpolar_medium
    aa_3 = ["F", "W", "Y"]  # aromatic
    aa_4 = ["S", "T", "C", "N", "Q"]  # polar
    aa_5 = ["K", "R", "H"]  # plus
    aa_6 = ["D", "E"]  # minus

    with open(f"{fasta_folder}\\{file}.fasta") as f:
        for line in f:
            if line[:1].strip() != ">":
                if line[number].strip() in aa_1:
                    return 1
                if line[number].strip() in aa_2:
                    return 2
                if line[number].strip() in aa_3:
                    return 3
                if line[number].strip() in aa_4:
                    return 4
                if line[number].strip() in aa_5:
                    return 5
                if line[number].strip() in aa_6:
                    return 6
                else:
                    return 0


def align_separate_sequences(fasta_file, pdb, destination):
    with open(f"{fasta_file}", "r") as file_input:
        with open(f"{destination}\\{pdb}.fasta", "w+") as file_output:
            start = f">{pdb}.pdb"
            end = ">"
            seq = ""
            for line in file_input:
                if line.strip() == start:
                    print(line.strip(), file=file_output)
                    for line2 in file_input:
                        if line2[:1].strip() == end:
                            break
                        seq = seq + line2.strip()
            print(seq.strip(), file=file_output)


def read_lenght(file, directory):
    with open(f"{directory}\\{file}.fasta") as f:
        for line in f:
            if line[:1] == ">":
                continue
            return len(line.strip())


def matrix_pdbs(list_of_pdbs):
    matrix = []
    for i in list_of_pdbs:
        for j in list_of_pdbs:
            if i != j or j != i:
                if f"{i}, {j}" in matrix or f"{j}, {i}" in matrix:
                    continue
                else:
                    matrix.append(f"{i}, {j}")
    return matrix


def check_if_done(project, structure_folder, list_struct):
    """simple prompt function to check if the structure folder contains everything the user wants"""
    print()
    pdb_present = set([file.replace(".pdb", "") for file in os.listdir(structure_folder) if
                       os.path.isfile(os.path.join(structure_folder, file))])

    useless_structures = [x for x in pdb_present if x not in list_struct]
    print(useless_structures)

    if set(list_struct + useless_structures) == set(pdb_present):
        print(f"All structures for the project {project} are present in {structure_folder}")
        return True
    else:
        for i in list_struct:
            if i in pdb_present:
                pass
            else:
                print(f"Missing structure : {i}.pdb")
        print(f"Please add missing structures to the folder {structure_folder}")
        continue_or_not("Are you finished?", "")
        check_if_done(project, structure_folder, list_struct)


def download_pdbs(liste, destination_folder):
    """download the files that needs to be downloaded, can ignore whose that does not exists (i.e. your structure)
    returns the list of structure that need to be placed manually"""
    structure_dl = []
    structure_to_dl = []

    for i in liste:
        if os.path.isfile(destination_folder + f"\\{i}.pdb"):
            print(f"File {i}.pdb already exists.")
            structure_dl.append(i)
        else:
            if is_pdb_valid(i):
                print(f"The pdb code {i} seems to be valid.")
                url_rcsb = f"https://files.rcsb.org/download/{i}.pdb"
                wget.download(url_rcsb, destination_folder)
                if os.path.isfile(destination_folder + f"\\{i}.pdb"):
                    print(f"File {i}.pdb downloaded successfully.")
                    structure_dl.append(i)
            else:
                print(f"The pdb code {i} doesn't seem to be valid, therefore it cannot be downloaded.")
                continue_or_not("Is it one of your personnal pdb files?",
                                "Cannot access this code on RCSB. Try to add it manually")
                structure_to_dl.append(i)

    if set(structure_dl + structure_to_dl) == set(liste):
        print(f"Structures downloaded successfully in the folder '{destination_folder}'")
        print()
    else:
        print("ERROR: Structures not downloaded successfully")
        exit(-1)

    return structure_to_dl


def continue_or_not(question, exit_message):
    userinp = input(f"{question} Y/N : ")
    if userinp.lower() == "y" or userinp.lower() == "yes":
        pass
    elif userinp.lower() == "n" or userinp.lower() == "no":
        print(f"{exit_message}")
        print("Exiting...")
        exit(-1)
    else:
        print("The answer must be Yes (Y) or No (N). Please retry")
        continue_or_not(question, exit_message)


def is_pdb_valid(pdb):
    """check if the pdb is valid by checking on the website rcsb"""
    autorized_characters = "abcdefghijklmnopqrstuvwxyz0123456789"
    ok_characters = 0
    urlt = f"https://www.rcsb.org/structure/{pdb}"
    resp = requests.get(urlt)
    if len(pdb) != 4:
        return False
    else:
        for each_character in pdb:
            if autorized_characters.count(each_character):
                ok_characters += 1
        if ok_characters == 4 and resp.status_code == 200:
            return True
        else:
            return False


def dictionary_coordinates_pdb(pdb_list, pdb_folder):
    pdb_dictionnary = {}
    for structure in pdb_list:
        with open(f"{pdb_folder}\\{structure}.pdb") as f:
            a_number = 0
            for line in f:
                if line[:4].strip() == "ATOM" and line[21].strip() == "A" and \
                        line[12:16].strip() == "CA" and line[16].strip() == "":

                    a_number += 1
                    pdb_dictionnary.update({f"{structure}, {a_number}": f"{line[30:38].strip()}, "
                                                                        f"{line[38:46].strip()}, "
                                                                        f"{line[46:54].strip()},"
                                                                        f"{line[22:26].strip()}"})

                elif line[:4].strip() == "ATOM" and line[21].strip() == "A" and \
                        line[12:16].strip() == "CA" and line[16].strip() == "A":

                    a_number += 1
                    pdb_dictionnary.update({f"{structure}, {a_number}": f"{line[30:38].strip()}, "
                                                                        f"{line[38:46].strip()}, "
                                                                        f"{line[46:54].strip()},"
                                                                        f"{line[22:26].strip()}"})
    return pdb_dictionnary


def how_many_pdb(directory):
    return len([name for name in os.listdir(directory) if os.path.isfile(os.path.join(directory, name))])


def get_list_from_ali(file):
    pdb_liste = []
    with open(file) as file_input:
        for line in file_input:
            if line[:1] == ">":
                pdb_liste.append(line[1:].strip().replace(".pdb", ""))

    return pdb_liste


def download_wget(adress, filename, destination):
    print(f"Downloading {filename} from {adress} into folder {destination}")
    wget.download(adress + f"{filename}", destination)
    if os.path.isfile(destination + f"\\{filename}"):
        print(f"File {filename} downloaded successfully")
    else:
        print(f"File {filename} not downloaded")
    return destination + f"\\{filename}"


def create_folder(path, folder):
    """create a named folder in a given path"""
    directory = path + f"\\{folder}"
    try:
        os.mkdir(directory)
    except OSError:
        print(f"'{folder}' folder exists already.")
    else:
        print(f"'{folder}' folder created sucessfully")
    return directory


def create_folder_destructive(path, folder):
    """create a named folder in a given path, if it allready exists, delete its files"""
    directory = path + f"\\{folder}"
    try:
        os.mkdir(directory)
    except OSError:
        print(f"'{folder}' folder exists already. Erasing its content.")
        delete_files_in_folder(directory)
    else:
        print(f"'{folder}' folder created sucessfully")
    return directory


def delete_files_in_folder(path):
    """delete all files in a given folder. be careful not to use anywhere,
    its like a biological bomb for humans but for files in a folder"""
    files = [name for name in os.listdir(path) if os.path.isfile(os.path.join(path, name))]
    for file in files:
        path_to_file = f"{path}" + f"\\{file}"
        os.remove(path_to_file)


# Main execution
if __name__ == '__main__':
    main()
