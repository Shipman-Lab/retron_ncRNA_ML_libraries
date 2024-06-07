import os
import shutil
import gzip

def find_filepaths(directory, substring):
    """
    finds all filepaths that contains a substring
    in a directory

    inputs:
    directory:		a directory (string)
    substring: 		the substring to find in the lower files

    output:
    filepaths:		all filepaths containing that substring in a dir (list of strings)
    """
    filepaths = []
    for root, dirs, files in os.walk(directory):
        if len(dirs) != 0:
            for subdir in dirs:
                subdir_files = os.listdir(os.path.join(root, subdir))
                for name in subdir_files:
                    if substring in name:
                        full_path = os.path.join(root, subdir, name)
                        filepaths.append(full_path)

        else:
            subdir_files = os.listdir(root)
            for name in subdir_files:
                if substring in name:
                    full_path = os.path.join(root, name)
                    filepaths.append(full_path)

    return filepaths

def unzip_if_necessary(full_file_paths, run_id):
    """
    checks if there is an unzipped copy of a run & unzips the gz if not

    inputs:
    full_file_paths: list of strings, each of which is a filepath

    output:
    unzipped_file_path: string, full filepath to unzipped file
    """
    # check that there are only 1 or 2 files in full filepath (aka
    # only a .gz & .fastq, not any duplicates or mismatches)
    if len(full_file_paths) > 2:
        raise ValueError("For run ID: %s has more than 2 matched files" % run_id)
    elif len(full_file_paths) == 0:
        raise ValueError("Run ID %s doesn't have a matching fastq file" % run_id)
    elif (len(full_file_paths) == 1) & ("fastq.gz" in full_file_paths[0]):
        with gzip.open(full_file_paths[0], 'rb') as f_in:
            with open(full_file_paths[0][:-3], 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
                unzipped_file_path = full_file_paths[0][:-3]
    elif len(full_file_paths) == 2:
        for path in full_file_paths:
            if (path[-5:] == "fastq") | (path[-5:] == "fasta"):
                unzipped_file_path = path
    elif (len(full_file_paths) == 1) & (".fasta" in full_file_paths[0]):
        unzipped_file_path = full_file_paths[0]
    else:
        raise ValueError("Unknown case in unzip_if_necessary")
    return unzipped_file_path