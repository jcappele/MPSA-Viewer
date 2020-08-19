# Import
import os
import zipfile

import requests
import wget

### Input data - Make sure everything in here is correct else the program will probably not work
# Name of the project
project_name = "Delta_boucle"
# List of the name of your structure files - Make sure to not use '.pdb' at the end
list_of_structure = ['2cc3_db', '3ub1_CD_db', '3wz3_db', '3wz4_db', '4akz_db', '4ec6_db', '4jf8_db', '4kz1_db',
                     '4lso_db', '4mei_db', '4nhf_db', '4o3v_db', '5aiw_db', '5cnl_db', '5i97_db', '6iqt_db', '6zgn_db']


### Main function
def main():
    """use of auxiliary functions to get the job done!
    in this case, this function will prepare the wanted structures
    for deposition to the mTM-align servers at yanglab.nankai.edu.cn"""
    print(f"Initialization of project named '{project_name}'")
    project_folder = create_folder(os.getcwd(), project_name)
    structure_folder = create_folder(os.getcwd(), "all_structures")
    msta_results = create_folder(project_folder, "MSTA_raw_results")
    processed_folder = create_folder(project_folder, "processed_pdbs")
    zip_folder = create_folder(project_folder, "zip_files")

    structure_left_behind = download_pdbs(list_of_structure, structure_folder)
    for i in structure_left_behind:
        if os.path.isfile(structure_folder + f"\\{i}.pdb"):
            pass
        else:
            print(f"The structure folder is lacking {i}.pdb, "
                  f"please add it in the folder {project_name}\\structure_folder")

    if check_if_done(project_name, structure_folder, list_of_structure):
        pass

    print()

    print("Transforming SeMethionine (MSE) to Methionines (MET). This step is essential for the alignement")
    mse_to_met(structure_folder, processed_folder, list_of_structure)
    print()
    print("All structures have been processed")
    print()
    print(f"Zipping all the processed structures into {project_name}\\{zip_folder}")
    create_zip(processed_folder, zip_folder, project_name, list_of_structure)
    print()
    print("The archive is ready to be uploaded to the mTM-server")
    print()
    print(f"Once the results are available, please place the fasta "
          f"alignement and the structure in the folder {msta_results} ")


### Auxiliary functions
def create_zip(f_folder, destination, project_n, structure_list):
    """create a zip file in a destination folder, based on files in an original folder
    the name of the archive is the one of the project + '.zip'"""
    zip_name = f"{destination}\\{project_n}.zip"
    with zipfile.ZipFile(zip_name, 'w') as zipf:
        for pdb in structure_list:
            print(f"Writing {pdb} to the archive {project_n}.zip")
            zipf.write(f"{f_folder}\\{pdb}.pdb", os.path.basename(f"{f_folder}\\{pdb}.pdb"))


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
                url = f"https://files.rcsb.org/download/{i}.pdb"
                wget.download(url, destination_folder)
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


def is_pdb_valid(pdb):
    """check if the pdb is valid by checking on the website rcsb"""
    autorized_characters = "abcdefghijklmnopqrstuvwxyz0123456789"
    ok_characters = 0
    url = f"https://www.rcsb.org/structure/{pdb}"
    resp = requests.get(url)
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


def delete_files_in_folder(folder):
    """delete all files in a given folder. be careful not to use anywhere,
    its like a biological bomb for humans but for files in a folder"""
    path = os.getcwd() + f"\\{folder}"
    files = [name for name in os.listdir(path) if os.path.isfile(os.path.join(path, name))]
    for file in files:
        path_to_file = os.getcwd() + f"\\{folder}" + f"\\{file}"
        os.remove(path_to_file)


### Main function execution
if __name__ == '__main__':
    main()
