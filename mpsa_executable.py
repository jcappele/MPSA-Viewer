#    MPSA Viewer, a simple and intuitive way to emphasize your multiple 
#    proteins alignements.

#    Copyright (C) 2020  Julien Cappèle, 
#                        Claude Didierjean, 
#                        Frédérique Favier
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.


import math
import os
import sys
import time
import subprocess
import zipfile
import urllib.request
import urllib.error
import shutil
import codecs

# Variables
version = "alpha 0.3"
config_file_name = "configuration_file"
parameters_file_name = "parameters"
structure_file_name = "structure_list"


# Main function
def main(begin: bool):
    if begin:
        program_startup()

    projects_folder = create_folder(os.getcwd(), "MPSA_Projects", False)

    user_input = ""
    while user_input != "exit":
        user_input = input("> ")

        # HELP
        if user_input == "help":
            get_help()
            main(False)

        # EXIT... again
        elif user_input == "exit":
            exit(0)

        # NEW PROJECT
        elif user_input.split(" ")[0] == "new_project":

            template_present = False
            template_name = ""
            try:
                project_name: str = user_input.split(" ")[1]
            except IndexError:
                project_name = ""
                print("You must specify a name to your new project.")
                main(False)

            try:
                template_name: str = user_input.split(" ")[2]
                template_present = True
            except IndexError:
                pass

            if os.path.exists(projects_folder + os.sep + project_name):
                print("This project already exists. Please choose another name")
                main(False)

            project_folder = create_folder(projects_folder, project_name, True)
            structure_folder = create_folder(project_folder, "Original_Structures", True)
            fastas_folder = create_folder(project_folder, "Alignement", True)
            pymol_session = create_folder(project_folder, "PyMOL_Session", True)
            pymol_script_folder = create_folder(project_folder, "PyMOL_Script", True)
            processed_folder = create_folder(project_folder, "Processed_Structures", True)
            aligned_structures_folder = create_folder(project_folder, "Aligned_Structures", True)
            zip_folder = create_folder(project_folder, "Zip_Files", True)
            msta_results_folder = create_folder(project_folder, "MSTA_result", True)

            if not template_present:
                generate_project_configuration(project_name, project_folder, structure_folder, fastas_folder,
                                               pymol_session, pymol_script_folder, processed_folder,
                                               aligned_structures_folder, msta_results_folder, zip_folder)
                generate_project_parameters(project_folder)

            if template_present:
                print(f"Using project {template_name} configuration files to generate project {project_name} files")
                try:
                    copy_config_parameters_from_template(project_name, template_name, projects_folder)
                except FileNotFoundError:
                    try:
                        shutil.rmtree(projects_folder + os.sep + project_name)
                    except FileNotFoundError:
                        print("ERROR2. This is getting out of hand. PANIC MODE")
                    print("ERROR. One of the configuration files is missing from the template folder.")
                    main(False)

            generate_structure_list(project_folder)

            time.sleep(1)
            print(f"Project {project_name} successfully generated. "
                  f"You may proceed to modify structure, configuration and parameters files")

        # EXECUTE ZIP CREATION
        elif user_input.split(" ")[0] == "zip":
            try:
                project_name: str = user_input.split(" ")[1]
            except IndexError:
                project_name = ""
                print("You must specify the project name.")
                main(False)

            try:
                config_dict = load_config_file(project_name, projects_folder, config_file_name)
            except FileNotFoundError:
                config_dict = {}
                print("ERROR: Project files not found. Please check the project name and/or folder files.")
                main(False)

            structure_folder = config_dict.get("structure_folder")
            processed_folder = config_dict.get("processed_folder")
            zip_folder = config_dict.get("zip_folder")

            list_of_structure_in_structure_file = load_structure_file(project_name, projects_folder)
            dict_of_structure = {}
            for structure in list_of_structure_in_structure_file:
                if structure[-2] == "_":
                    dict_of_structure.update({structure[:-2]: structure[-1]})
                else:
                    dict_of_structure.update({structure: "A"})

            list_of_structures_to_dl = [structure for structure in dict_of_structure.keys()]

            structure_left_behind = download_pdbs(list_of_structures_to_dl, structure_folder)

            for i in structure_left_behind:

                if os.path.isfile(structure_folder + os.sep + f"{i}.pdb"):
                    pass
                else:
                    print(f"The structure folder is lacking {i}.pdb, "
                          f"please add it in the folder {project_name}{os.sep}structure_folder")

            if check_if_done(project_name, structure_folder, list_of_structures_to_dl):
                pass

            for file in os.listdir(structure_folder):
                mse_to_met(structure_folder, processed_folder, file, dict_of_structure.get(file[:-4]))

            create_zip(processed_folder, zip_folder, project_name, list_of_structures_to_dl)
            print("Zip file created. You can submit it to mTM-align server")

        # DOWNLOAD MSTA_results
        elif user_input.split(" ")[0] == "fetch":
            try:
                project_name: str = user_input.split(" ")[1]
            except IndexError:
                project_name = ""
                print("You must specify the project name.")
                main(False)

            try:
                url_download: str = user_input.split(" ")[2]
            except IndexError:
                url_download = ""
                print("You must specify the download url.")
                main(False)

            try:
                config_dict = load_config_file(project_name, projects_folder, config_file_name)
            except FileNotFoundError:
                config_dict = {}
                print("ERROR: Project files not found. Please check the project name and/or folder files.")
                main(False)

            msta_results_folder = config_dict.get("msta_results_folder")

            print("""    WARNING: for whatever reason, downloading results
    can take several minutes. You might want to 
    download them manually on your favorite
    webbrowser. More details in the documentation.""")
            print()
            print("Downloading files... Please wait.")
            urllib.request.urlretrieve(url_download + "/new.pdb",
                                       msta_results_folder + os.sep + "structure_alignement.pdb")
            urllib.request.urlretrieve(url_download + "/seq.fasta",
                                       msta_results_folder + os.sep + "seq_alignement.fasta")
            print("Download completed.")

        # EXECUTE MSTA ANALYSIS
        elif user_input.split(" ")[0] == "analysis":
            try:
                project_name: str = user_input.split(" ")[1]
            except IndexError:
                project_name = ""
                print("You must specify the project name.")
                main(False)

            try:
                parameters_dict = load_config_file(project_name, projects_folder, parameters_file_name)
            except FileNotFoundError:
                parameters_dict = {}
                print("ERROR: Project files not found. Please check the project name and/or folder files.")
                main(False)

            try:
                config_dict = load_config_file(project_name, projects_folder, config_file_name)
            except FileNotFoundError:
                config_dict = {}
                print("ERROR: Project files not found. Please check the project name and/or folder files.")
                main(False)

            msta_results_folder = config_dict.get("msta_results_folder")
            processed_folder = config_dict.get("processed_folder")
            fasta_folder = config_dict.get("fasta_folder")
            pymol_session = config_dict.get("pymol_session")
            pymol_script_folder = config_dict.get("pymol_script_folder")
            aligned_structures_folder = config_dict.get("aligned_structures_folder")
            pymol_path = config_dict.get("pymol_path").replace(r"\"", "/")

            time_multiplier = parameters_dict.get("time_multiplier")
            coloration_method = parameters_dict.get("coloration_method")
            coloration_cutoff = parameters_dict.get("coloration_cutoff")
            background_color = rgb_to_hexa_str(parameters_dict.get("background_color"))
            default_coloration_backbone = bool(parameters_dict.get("default_coloration_backbone"))
            initial_color_backbone = rgb_to_hexa_str(parameters_dict.get("initial_color_backbone"))
            final_color_backbone = rgb_to_hexa_str(parameters_dict.get("final_color_backbone"))
            number_of_variantes = parameters_dict.get("number_of_variantes")
            not_aligned_backbone = parameters_dict.get("not_aligned_backbone")
            aminoacids_representation = parameters_dict.get("aminoacids_representation")
            default_coloration_aminoacids = bool(parameters_dict.get("default_coloration_aminoacids"))
            identical_aminoacids_color = rgb_to_hexa_str(parameters_dict.get("identical_aminoacids_color"))
            similar_aminoacids_color = rgb_to_hexa_str(parameters_dict.get("similar_aminoacids_color"))
            default_representation_aminoacids = bool(parameters_dict.get("default_representation_aminoacids"))
            identical_aminoacids_representation = parameters_dict.get("identical_aminoacids_representation")
            similar_aminoacids_representation = parameters_dict.get("similar_aminoacids_representation")
            ligand_representation = parameters_dict.get("ligand_representation")

            if default_coloration_backbone:
                initial_color_backbone = "10, 10, 255"
                final_color_backbone = "255, 10, 10"
                number_of_variantes = 10
                not_aligned_backbone = rgb_to_hexa_str("50, 50, 50")

            value_color_dict = color_madness(initial_color_backbone, final_color_backbone, number_of_variantes,
                                             coloration_cutoff)

            if default_coloration_aminoacids:
                identical_aminoacids_color = rgb_to_hexa_str("0, 255, 0")
                similar_aminoacids_color = rgb_to_hexa_str("10, 10, 10")

            if default_representation_aminoacids:
                identical_aminoacids_representation = "sticks"
                similar_aminoacids_representation = "lines"
                ligand_representation = "sticks"

            list_dir = os.listdir(msta_results_folder)
            if len(list_dir) != 2:
                print("ERROR: There are more than 2 files in the result folder. Please remove unnecessary ones.")
                main(False)

            pdb_file_present = False
            seq_file_present = False
            name_fasta_file = ""
            name_pdb_file = ""

            for item in list_dir:
                if item[-5:] == "fasta":
                    seq_file_present = True
                    name_fasta_file = item.strip()
                if item[-3:] == "pdb":
                    pdb_file_present = True
                    name_pdb_file = item.strip()
            if seq_file_present and pdb_file_present:
                print("Results OK!")
            else:
                print("ERROR: One of the result file is not in a compatible format. Please check result folder")
                main(False)

            path_to_fasta = (msta_results_folder + os.sep + name_fasta_file)
            path_to_pdb = (msta_results_folder + os.sep + name_pdb_file)
            pdb_list = get_list_from_ali(path_to_fasta)

            if coloration_method == "RMSD" and len(pdb_list) < 3:
                print()
                print("WARNING. RMSD method cannot be used with less than 3 structures.\n"
                      "         Changing coloration method to highest distances")
                print()
                time.sleep(2)
                coloration_method = "HD"

            message_01 = "Aligning those structures: "
            for pdb in pdb_list:
                message_01 += pdb + " "

            with open(pymol_script_folder + os.sep + f"super{project_name}.pml", "w") as output:
                print(f"""reinitialize
                        load {path_to_pdb}""", file=output)
                for x in pdb_list:
                    print(f"""load {processed_folder}{os.sep}{x}.pdb
                    super {x}, {name_pdb_file.replace(".pdb", "")}""", file=output)
                for x in pdb_list:
                    print(f"""save {aligned_structures_folder}{os.sep}{x}.pdb, {x}""", file=output)
                print("""quit""", file=output)
            print(f"Waiting for {len(pdb_list) * 1.0 * float(time_multiplier)} "
                  f"seconds, time to superimpose all structures with PyMOL.")

            try:
                subprocess.call([pymol_path, pymol_script_folder + os.sep + f"super{project_name}.pml"])
            except FileNotFoundError or PermissionError:
                print("PyMOL path not found or incorrect. "
                      "Make sure you added the correct path in the configuration file.")
                main(False)

            time.sleep(len(pdb_list) * 1.0 * float(time_multiplier))

            aligned_files = [file[:4] for file in os.listdir(aligned_structures_folder)]
            if aligned_files != pdb_list:
                print("ERROR. It seems that PyMOL didn't complete his alignement in time.\n"
                      "Make sure to increase the time_multipler value in the parameters file.")
                main(False)
            else:
                print("Alignement complete!")

            pdb_dict = dictionary_coordinates_pdb(pdb_list, aligned_structures_folder)

            for i in pdb_list:
                print(f"Creation of {i}.fasta")
                align_separate_sequences(path_to_fasta, i, fasta_folder)
            hm_fasta = how_many_pdb(fasta_folder)
            print("Dictionary of structures and coordinates generated correctly")

            matrix_pdb = matrix_pdbs(pdb_list)
            print("Matrices of structures generated correctly")
            align_len = read_lenght(pdb_list[0], fasta_folder)

            with open(pymol_script_folder + os.sep + f"{project_name}.pml", "w") as f:
                print("""reinitialize""", file=f)
                for n in pdb_list:
                    print(f"""load {aligned_structures_folder}{os.sep}{n}.pdb""", file=f)
                print(f"""color {not_aligned_backbone}, all
                    select water, resn HOH
                    remove water""", file=f)

            with open(pymol_script_folder + os.sep + f"{project_name}_aminoacids.pml", "w"):
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
                    if check_pos_ifchar(chara, structure, fasta_folder):
                        number += 1
                    aa_type.update(give_char(chara, structure, fasta_folder))
                    aa_class.update(str(check_pos_whichclass(chara, structure, fasta_folder)))

                if number >= hm_fasta:
                    dict_dist = {}
                    for tuplex in range(len(matrix_pdb)):
                        x = matrix_pdb[tuplex]
                        une, deux = x.split(",")

                        a = une.strip()
                        b = deux.strip()

                        # structure 1 - a
                        seq_a = fetch_sequence_in_fasta(a, fasta_folder)
                        dash_a = count_dash(seq_a, 0, position)
                        calculus_a = position - dash_a
                        x1, y1, z1, seqres1 = xyzseq_pdb(pdb_dict, a, calculus_a)

                        # structure 2 - b
                        seq_b = fetch_sequence_in_fasta(b, fasta_folder)
                        dash_b = count_dash(seq_b, 0, position)
                        calculus_b = position - dash_b
                        x2, y2, z2, seqres2 = xyzseq_pdb(pdb_dict, b, calculus_b)

                        dict_dist.update(
                            {f"{a}, {seqres1}, {b}, {seqres2}": f"{position}, {distance(x1, y1, z1, x2, y2, z2)}"})

                    plus_grande_valeur.append(max(dict_dist.values()))
                    result, pos = 0, 0

                    if coloration_method == "HD":
                        value = max(dict_dist.values())
                        pos, result = value.split(",")
                        pos, result = int(pos.strip()), float(result.strip())

                    elif coloration_method == "RMSD":
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
                    with open(pymol_script_folder + os.sep + f"{project_name}.pml", "a") as f:
                        with open(pymol_script_folder + os.sep + f"{project_name}_aminoacids.pml", "a") as f2:
                            for i in pdb_list:
                                for item in value_color_dict:
                                    if result <= value_color_dict.get(item):
                                        print_in_script(i, pos, item, f, pdb_dict, fasta_folder)
                                        break
                                    else:
                                        pass

                            for i in pdb_list:
                                if len(aa_class) == 1 and len(aa_type) == 1:
                                    print_in_script_aminoacids(i, position, identical_aminoacids_color, f2, pdb_dict,
                                                               fasta_folder, identical_aminoacids_representation)
                                if len(aa_class) == 1 and len(aa_type) > 1:
                                    print_in_script_aminoacids(i, position, similar_aminoacids_color, f2, pdb_dict,
                                                               fasta_folder, similar_aminoacids_representation)

                else:
                    pass

            with open(pymol_script_folder + os.sep + f"{project_name}.pml", "a") as f:
                if aminoacids_representation:
                    print(f"""run {pymol_script_folder}{os.sep}{project_name}_aminoacids.pml""", file=f)
                print(f"""set cartoon_flat_sheets, 1
set cartoon_side_chain_helper, 1
set cartoon_fancy_helices, 1
select Identical_Residues, color {identical_aminoacids_color}
hide (Identical_Residues)
show {identical_aminoacids_representation}, Identical_Residues
select SameClass_Residues, color {similar_aminoacids_color}
hide (SameClass_Residues)
show {similar_aminoacids_representation}, SameClass_Residues
hide (hetatm)
show {ligand_representation}, hetatm
orient all
bg_color {background_color}
util.cnc
deselect
save {pymol_session}{os.sep}{project_name}.pse""", file=f)

            print()
            print("Opening PyMOL session")

            subprocess.call([pymol_path, pymol_script_folder + os.sep + f"{project_name}.pml"])

        # HELP COMMAND
        elif user_input == "help exit":
            print(exit_pg.__doc__)

        # NO VIABLE COMMAND ENTERED
        else:
            print(f"Unknown command : {user_input}")


# Auxilliary commands
def get_help():
    print("""Here is the list of commands:
    
exit                                cancel the execution of the program


new_project <name> (<template>)     generate a new project folder
                                    <name> = name of your project
                                    (<template>) name of the template
                                    
                IF you specify a template name
                your new project will have the
                same configuration and parameters
                files. Not mandatory
                                
                                
zip <name>                          prepare the pdb files for deposition
                                    <name> = name of your project


fetch_results <name> <url>          download the MSTA results of a project
                                    <name> = name of your project
                                    <url> = your url to the completed 
                                    mTM-align result
                                
            WARNING: for whatever reason, downloading results
            can take several minutes. You might want to 
            download them manually on your favorite
            webbrowser. More details in the documentation.
                                
                                
analysis <name>                     execute the alignement analysis
                                    and create a pymol session based
                                    on your parameters                                    
                                    
                                    """)


def copy_config_parameters_from_template(project_name, template_name, projects_folder):
    with codecs.open(projects_folder + os.sep + template_name + os.sep + config_file_name + ".txt",
                     "r", encoding='utf-8') as template_config:
        with codecs.open(projects_folder + os.sep + project_name + os.sep + config_file_name + ".txt",
                         "w", encoding='utf-8') as project_config:
            for line in template_config:
                print(line, file=project_config)

    with codecs.open(projects_folder + os.sep + template_name + os.sep + parameters_file_name + ".txt",
                     "r", encoding='utf-8') as template_param:
        with codecs.open(projects_folder + os.sep + project_name + os.sep + parameters_file_name + ".txt",
                         "w", encoding='utf-8') as project_param:
            for line in template_param:
                print(line, file=project_param)


def mse_to_met(location_folder, destination_folder, file, chain):
    """change the selenomethionine in methionine for the pdb file
    write new file into a destination folder"""
    with open(f"{location_folder}{os.sep}{file}", "r") as file_input:
        print(f"Processing {file}")
        with open(f"{destination_folder}{os.sep}{file}", "w") as file_output:
            for line in file_input:
                if line[:6].strip() == "HETATM" or line[:4].strip() == "ATOM":
                    if line[21:22].strip() != chain:
                        continue
                if line[:6].strip() == "HETATM" and line[17:20].strip() == "MSE":
                    print(line.strip().replace("HETATM", "ATOM  ").replace("MSE", "MET").replace("SE", " S"),
                          file=file_output)
                else:
                    print(line.strip(), file=file_output)
        print(f"{file} successfully processed")


def load_config_file(project_name, projects_folder, file_name):
    configfile_dict = {}
    with codecs.open(projects_folder + os.sep + project_name + os.sep + f"{file_name}.txt",
                     encoding='utf-8') as config_file:
        lines = filter(None, (line.rstrip() for line in config_file))
        for line in lines:
            if not line.startswith("#"):
                try:
                    key, value = line.strip().split(" = ")
                except ValueError:
                    print(
                        f"ERROR. One or more values in the {file_name}.txt file is empty/uncorrectly filled. "
                        f"Please check it and retry.")
                    main(False)

                configfile_dict.update({key: value})

    return configfile_dict


def load_structure_file(project_name, projects_folder):
    structure_list = []
    with codecs.open(projects_folder + os.sep + project_name + os.sep + f"structure_list.txt",
                     encoding='utf-8') as file:
        lines = filter(None, (line.rstrip() for line in file))
        for line in lines:
            if not line.startswith("#"):
                structure_list.append(line.strip())

    return structure_list


def get_projects(path_to_projects_folder):
    projects = []
    for root, dirs, files in os.walk(path_to_projects_folder):
        for file in files:
            if file.endswith(f"{config_file_name}.txt"):
                projects.append(file.replace(f"{config_file_name}.txt", ""))
    if len(projects) == 0:
        print("No project detected.")
    else:
        string = "These projects are available : "
        for project in projects:
            string += (project + " ")
        print(string)


def generate_structure_list(project_folder):
    with open(project_folder + os.sep + f"{structure_file_name}.txt", "w") as structure_list:
        print("""# List of structures. One structure per line. No special characters nor spaces in pdb names.
# You can specify the chain at the end of the line like this : 3ec6_B or mypdbfile_B""", file=structure_list)


def generate_project_configuration(project_name, project_folder, structure_folder, fastas_folder, pymol_session,
                                   pymol_script_folder, processed_folder, aligned_structures_folder,
                                   msta_results_folder, zip_folder):
    with open(project_folder + os.sep + f"{config_file_name}.txt", "w") as configuration_file:
        print(f"""# CONFIGURATION FILE
# PyMOL path ex: C:\\Users\\User\\PyMOL\\PyMOLWin.exe 
# Make sure the separators are compatible with your operating system
pymol_path = 


#  GENERATED BY THE PROGRAM, DO NOT TOUCH UNLESS YOU KNOW WHAT YOU ARE DOING
project_name = {project_name}
project_folder = {project_folder}
msta_results_folder = {msta_results_folder}
structure_folder = {structure_folder}
zip_folder = {zip_folder}
fasta_folder = {fastas_folder}
pymol_session = {pymol_session}
pymol_script_folder = {pymol_script_folder}
processed_folder = {processed_folder}
aligned_structures_folder = {aligned_structures_folder}""", file=configuration_file)


def generate_project_parameters(project_folder):
    with open(project_folder + os.sep + f"{parameters_file_name}.txt", "w") as parameters_list:
        print("""#  PARAMETERS FOR MPSA VIEWER

#  PyMOL alignement time multiplier. 
#  If you encounter an error during PyMOL aligning the structures, try increasing this value
#  Default is 1
time_multiplier = 1
        
#  Coloration method to use for the superposition. Valid values : "RMSD" or "HD"
#  WARNING : RMSD method is only usable with 3 or more structures
coloration_method = RMSD

#  Coloration distance cut-off. Default value for RMSD = 2, default value for HD = 10
coloration_cutoff = 2

#  Coloration of the background in PyMOL. Default is 255, 255, 255
background_color = 255, 255, 255

#  Use default coloration for the backbone ? Valid values : True or False
default_coloration_backbone = True

#  Specify your desired colors for the backbone in RGB format, and the number of variantes you want (max = 100)
#  This will be ignored if default_coloration_backbone is set to True
#  For the best visual representation, try to set one color below 50 in both initial and final colors. 
#  You could end up with grey shades if you don't. Exemple of worst case scenario : 255, 125, 0 and 0, 125, 255
initial_color_backbone = 10, 10, 255
final_color_backbone = 255, 10, 10
number_of_variantes = 10
not_aligned_backbone = 50, 50, 50

#  Aminoacid representation in the PyMOL session. Valid values : True or False
aminoacids_representation = True

#  Use default colors for aminoacids ? Valid values : True or False
#  This will be ignored if aminoacids_representation is set to False
default_coloration_aminoacids = True

#  Specify your desired colors for aminoacids in RGB format.
#  This will be ignored if default_coloration_aminoacids is set to True
identical_aminoacids_color = 0, 255, 0
similar_aminoacids_color = 10, 10, 10

#  Use default representation for aminoacids and ligand ? Valid values : True or False
default_representation_aminoacids = True

#  Specify your disired representation for aminoacids and ligand. Valid values : sticks, lines, spheres, dots
#  This will be ignored if default_representation_aminoacids is set to True
identical_aminoacids_representation = sticks
identical_aminoacids_representation = lines
ligand_representation = sticks
""", file=parameters_list)


def exit_pg():
    """This command will exit the program"""
    sys.exit()


def program_startup():
    print(f"""
=======================================================================
============================= MPSA_VIEWER =============================
=======================================================================
                          version {version}
                          
Copyright (C) 2020 Julien Cappèle, Claude Didierjean, Frédérique Favier
This program comes with ABSOLUTELY NO WARRANTY.
This is free software, and you are welcome to redistribute it
under certain conditions.

Type "help" to get the list of commands""")


def create_zip(f_folder, destination, project_n, structure_list):
    """create a zip file in a destination folder, based on files in an original folder
    the name of the archive is the one of the project + '.zip'"""
    zip_name = f"{destination}{os.sep}{project_n}.zip"
    with zipfile.ZipFile(zip_name, 'w') as zipf:
        for pdb in structure_list:
            print(f"Writing {pdb} to the archive {project_n}.zip")
            zipf.write(f"{f_folder}{os.sep}{pdb}.pdb", os.path.basename(f"{f_folder}{os.sep}{pdb}.pdb"))


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


def check_if_done(project, structure_folder, list_struct):
    """simple prompt function to check if the structure folder contains everything the user wants"""
    list_struct = [elem[:-2] for elem in list_struct if elem[-2] == "_"]
    pdb_present = set([file.replace(".pdb", "") for file in os.listdir(structure_folder) if
                       os.path.isfile(os.path.join(structure_folder, file))])

    useless_structures = [x for x in pdb_present if x not in list_struct]

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
        print(i)
        print(destination_folder + os.sep + f"{i}.pdb")
        print(os.path.isfile(destination_folder + os.sep + f"{i}.pdb"))
        if os.path.isfile(destination_folder + os.sep + f"{i}.pdb"):
            print(f"File {i}.pdb already exists.")
            structure_dl.append(i)
        else:
            if is_pdb_valid(i):
                print(f"The pdb code {i} seems to be valid.")
                url = f"https://files.rcsb.org/download/{i}.pdb"
                print(f"Downloading structure {i}")
                urllib.request.urlretrieve(url, destination_folder + os.sep + f"{i}.pdb")
                print(f"Download successful")
                if os.path.isfile(destination_folder + os.sep + f"{i}.pdb"):
                    print(f"File {i}.pdb downloaded successfully.")
                    structure_dl.append(i)
            else:
                print(f"The pdb code {i} doesn't seem to be valid, therefore it cannot be downloaded.")
                continue_or_not("Is it one of your personnal pdb files?",
                                "Cannot access this code on RCSB. Try to add it manually")
                structure_to_dl.append(i)

    if set(structure_dl + structure_to_dl) == set(liste):
        print(f"Structures downloaded successfully in the folder '{destination_folder}'")
    else:
        print("ERROR: Structures not downloaded successfully")
        exit(-1)

    return structure_to_dl


def is_pdb_valid(pdb):
    """check if the pdb is valid by checking on the website rcsb"""
    if pdb[-2] == "_":
        url = f"https://www.rcsb.org/structure/{pdb[:-2]}"
    else:
        url = f"https://www.rcsb.org/structure/{pdb}"

    try:
        status_code = urllib.request.urlopen(url).getcode()
    except urllib.error:
        status_code = 404

    if status_code == 200:
        return True
    else:
        return False


def create_folder(path_of_folder, folder, verbose: bool):
    """create a named folder in a given path"""
    directory = path_of_folder + os.sep + f"{folder}"
    try:
        os.mkdir(directory)
    except OSError:
        if verbose:
            print(f"'{folder}' folder exists already.")
    else:
        if verbose:
            print(f"'{folder}' folder created sucessfully")
    return directory


def delete_files_in_folder(folder):
    """delete all files in a given folder. be careful not to use anywhere,
    its like a biological bomb for humans but for files in a folder"""
    path_of_flder = os.getcwd() + os.sep + f"{folder}"
    files = [name for name in os.listdir(path_of_flder) if os.path.isfile(os.path.join(path_of_flder, name))]
    for file in files:
        path_to_file = os.getcwd() + os.sep + f"{folder}" + os.sep + f"{file}"
        os.remove(path_to_file)


def color_madness(color_initial, color_final, how_many_variantes, cutoff):
    initial_red, initial_green, initial_blue = str(color_initial).split(",")
    final_red, final_green, final_blue = str(color_final).split(",")
    red_dif, green_dif, blue_dif = float(initial_red) - float(final_red), float(initial_green) - float(
        final_green), float(initial_blue) - float(final_blue)

    variantes_rgb = []
    for i in range(how_many_variantes):
        variante_red = int(float(initial_red) - (float(red_dif) * (i / (how_many_variantes - 1))))
        variante_green = int(float(initial_green) - (float(green_dif) * (i / (how_many_variantes - 1))))
        variante_blue = int(float(initial_blue) - (float(blue_dif) * (i / (how_many_variantes - 1))))
        variantes_rgb.append([variante_red, variante_green, variante_blue])

    variante_hexadecimal = [rgb_to_hexa_list(variante) for variante in variantes_rgb]

    cutoff_one_size = int(cutoff) / int(how_many_variantes)
    dict_variantes = {}
    cutoff_state = 1

    for item in variante_hexadecimal:
        dict_variantes.update({item: cutoff_one_size * cutoff_state})
        cutoff_state += 1

    return dict_variantes


def rgb_to_hexa_list(rgb: list):
    return f"0x{rgb[0]:02x}{rgb[1]:02x}{rgb[2]:02x}"


def rgb_to_hexa_str(rgb: str):
    rgb2 = [int(item.strip()) for item in rgb.split(",")]
    return f"0x{rgb2[0]:02x}{rgb2[1]:02x}{rgb2[2]:02x}"


def rmsd_list(list_val):
    size = len(list_val)
    mean_intermediate = 0
    for i in list_val:
        mean_intermediate += float(i)
    mean = mean_intermediate / size

    calculus = 0
    for i in list_val:
        calculus += ((float(i) - mean) ** 2)
    if size > 1:
        the_result = math.sqrt(calculus / (size - 1))
    else:
        the_result = math.sqrt(calculus / size)

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
    print(f"""show {representation}, {i} and resi {seqres}
color {color}, {i} and resi {seqres} and not (name N+O+C+CA)""", file=file)


def print_in_script(pdb, position, color, file, pdb_dictionnary, fasta_folder):
    seq = fetch_sequence_in_fasta(pdb, fasta_folder)
    dash = count_dash(seq, 0, position)
    true_pos = position - dash
    x, y, z, seqres = xyzseq_pdb(pdb_dictionnary, pdb, true_pos)
    i = pdb.replace(".pdb", "")
    print(f"""color {color}, {i} and resi {seqres} and (name c+o+n+ca)""", file=file)


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


def xyzseq_pdb(dic, pdb, number):
    x, y, z, nseq = (dic.get(f"{pdb}, {number}")).split(',')
    return float(x), float(y), float(z), int(nseq)


def distance(x1, y1, z1, x2, y2, z2):
    return math.sqrt(((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2))


def count_dash(string, start, end):
    return string.count("-", int(start), int(end))


def fetch_sequence_in_fasta(file, fasta_folder):
    with open(f"{fasta_folder}{os.sep}{file}.fasta") as f:
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
    with open(f"{fasta_folder}{os.sep}{file}.fasta") as f:
        for line in f:
            if line[:1].strip() != ">":
                if line[number].strip() in aminoacids:
                    return True
                else:
                    return False


def give_char(number, file, fasta_folder):
    with open(f"{fasta_folder}{os.sep}{file}.fasta") as f:
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

    with open(f"{fasta_folder}{os.sep}{file}.fasta") as f:
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
        with open(f"{destination}{os.sep}{pdb}.fasta", "w+") as file_output:
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
    with open(f"{directory}{os.sep}{file}.fasta") as f:
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


def dictionary_coordinates_pdb(pdb_list, pdb_folder):
    pdb_dictionnary = {}
    print(pdb_dictionnary)
    for structure in pdb_list:
        print(structure)
        with open(f"{pdb_folder}{os.sep}{structure}.pdb") as f:
            a_number = 0
            for line in f:
                if line[:4].strip() == "ATOM" and line[12:16].strip() == "CA" and line[16].strip() == "":

                    a_number += 1
                    pdb_dictionnary.update({f"{structure}, {a_number}": f"{line[30:38].strip()}, "
                                                                        f"{line[38:46].strip()}, "
                                                                        f"{line[46:54].strip()},"
                                                                        f"{line[22:26].strip()}"})

                elif line[:4].strip() == "ATOM" and line[12:16].strip() == "CA" and line[16].strip() == "A":

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


def create_folder_destructive(path_dir, folder):
    """create a named folder in a given path, if it allready exists, delete its files"""
    directory = path_dir + os.sep + f"{folder}"
    try:
        os.mkdir(directory)
    except OSError:
        print(f"'{folder}' folder exists already. Erasing its content.")
        delete_files_in_folder(directory)
    else:
        print(f"'{folder}' folder created sucessfully")
    return directory


if __name__ == '__main__':
    main(True)
