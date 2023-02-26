#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 25 09:35:32 2023

@author: os18o068
"""

import os
import csv
import pandas as pd
from datetime import datetime


def check_type(s):
    try:
        float(s)
        if "." in s:
            return float(s)
        else:
            return int(s)

    except ValueError:
        return s


def compare_dataframes(df1, df2):
    if len(df1.columns) != len(df2.columns):
        return False
    if not all(df1.columns == df2.columns):
        return False
    return True


def compare_data(dat1, dat2, met1, met2, which="both"):
    if which == "both":
        if compare_dataframes(dat1, dat2) and met1 == met2:
            return True

    elif which == "meta":
        raise NotImplementedError("The option which=meta is not available yet.")

    elif which == "data":
        raise NotImplementedError("The option which=both is not available yet.")

    else:
        return False


def read_csv(path, specs={}, comment="#"):
    # Read data
    with open(path) as fd:
        df = pd.read_csv(fd, delimiter=",", comment=comment)

    # Read in the comments and store them as a dictionary
    with open(path, "r") as f:
        comments = {}
        for line in f:
            if line.startswith(comment):
                key, value = line.strip(comment).strip().split(":")
                comments[key.strip()] = value.strip()

    metadata = {}
    current_key = None

    with open(path, "r") as f:
        for line in f:
            if line.startswith(comment):
                split_string = line.strip()[1:].strip().split(":")
                if line[1] != "!":
                    # check if line is key or key-value pair
                    # is key-value pair

                    if len(split_string) > 1 and split_string[1].strip() != "":
                        key_value = line.strip()[1:].strip().split(":")
                        if "--" in key_value[1]:
                            val = [
                                check_type(x.strip()) for x in key_value[1].split("--")
                            ]
                        else:
                            val = check_type(key_value[1].strip())

                        metadata[key_value[0].strip()] = val

                    # is not key-value pair
                    else:
                        current_key = split_string[0]
                        metadata[current_key] = {}

                else:
                    if current_key is None:
                        continue
                    key_value = line.strip()[2:].strip().split(":")
                    if "--" in key_value[1]:
                        val = [check_type(x.strip()) for x in key_value[1].split("--")]
                    else:
                        val = check_type(key_value[1].strip())

                    metadata[current_key][key_value[0].strip()] = val

    return df, metadata


def write_to_csv(
    data, target_file, conflict="ommit", meta={}, comment="#", check_compatibility=True
):

    if len(list(meta.keys())) != 0:
        add_meta = True

    else:
        add_meta = False

    # Check if the target file already exists
    file_exists = os.path.isfile(target_file)

    if file_exists:

        # If conflict is "ommit", create a buffer file with a time stamp
        if conflict == "ommit":
            timestamp = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
            buffer_file = f"{target_file.split('.')[0]}_buffer_{timestamp}.csv"
            new_df = pd.DataFrame(data)

            if add_meta:
                # Open a file for writing and write metadata as comments
                with open(buffer_file, "w", newline="") as csvfile:
                    writer = csv.writer(csvfile, delimiter=",")
                    for key, value in meta.items():
                        if isinstance(value, dict):
                            string = f"{comment} {key}:"
                            writer.writerow([string])
                            for k, v in value.items():
                                # check type of value
                                if isinstance(v, list):
                                    var_string = "--".join(str(x) for x in v)
                                else:
                                    var_string = str(v)
                                string = f"{comment}!{k}: {var_string}"
                                writer.writerow([string])
                        else:

                            if isinstance(value, list):
                                var_string = "--".join(str(x) for x in value)
                            else:
                                var_string = str(value)
                            string = f"{comment} {key}: {var_string}"
                            writer.writerow([string])

                    # Write DataFrame to the file without adding an empty row
                    new_df.to_csv(csvfile, index=False)

            else:
                new_df.to_csv(buffer_file, index=False)

            print(
                f"Warning: Target file {target_file} cannot be overwritten. Data has been saved to buffer file {buffer_file} instead. To merge buffer files use the merge_buffer_files() method. Only use this method if you know what you're doing!"
            )
            return

        # If file exists and conflict is "add", append new data to existing file
        if conflict == "add":
            # Read in existing data
            existing_df, existing_meta = read_csv(target_file)
            # chack compatibility between old and new data
            print("existing meta =", existing_meta)
            print("meta =", meta)
            compatible = compare_data(existing_df, data, existing_meta, meta)

            if compatible:
                combined_df = existing_df.append(data, ignore_index=True)
            else:
                if check_compatibility:
                    print(
                        f"Warning: Passed data format is inconsistent with original file. Pass check_compatibility = False to ignore this. Only use this option if you know what you're doing!"
                    )
                    timestamp = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
                    buffer_file = f"{target_file.split('.')[0]}_buffer_{timestamp}.csv"
                    new_df = pd.DataFrame(data)
                    new_df.to_csv(buffer_file, index=False)
                    print(
                        f"Warning: Target file {target_file} cannot be overwritten. Data has been saved to buffer file {buffer_file} instead. To merge buffer files use the merge_buffer_files() method."
                    )
                    return
                else:
                    combined_df = existing_df.append(data, ignore_index=True)

            data.to_csv(
                target_file, mode="a", index=False, header=not existing_df.index.any()
            )

        # If file exists and conflict is "overwrite", overwrite existing file with new data
        elif conflict == "overwrite":
            new_df = pd.DataFrame(data)

            if add_meta:
                # Open a file for writing and write metadata as comments
                with open(target_file, "w", newline="") as csvfile:
                    writer = csv.writer(csvfile, delimiter=",")
                    for key, value in meta.items():
                        if isinstance(value, dict):
                            string = f"{comment} {key}:"
                            writer.writerow([string])
                            for k, v in value.items():
                                # check type of value
                                if isinstance(v, list):
                                    var_string = "--".join(str(x) for x in v)
                                else:
                                    var_string = str(v)
                                string = f"{comment}!{k}: {var_string}"
                                writer.writerow([string])
                        else:

                            if isinstance(value, list):
                                var_string = "--".join(str(x) for x in value)
                            else:
                                var_string = str(value)
                            string = f"{comment} {key}: {var_string}"
                            writer.writerow([string])

                    # Write DataFrame to the file without adding an empty row
                    new_df.to_csv(csvfile, index=False)

            else:
                new_df.to_csv(target_file, index=False)

    # If file doesn't exist, create new file with new data
    else:
        new_df = pd.DataFrame(data)
        if add_meta:
            # Open a file for writing and write metadata as comments
            with open(target_file, "w", newline="") as csvfile:
                writer = csv.writer(csvfile, delimiter=",")
                for key, value in meta.items():
                    if isinstance(value, dict):
                        string = f"{comment} {key}:"
                        writer.writerow([string])
                        for k, v in value.items():
                            # check type of value
                            if isinstance(v, list):
                                var_string = "--".join(str(x) for x in v)
                            else:
                                var_string = str(v)
                            string = f"{comment}!{k}: {var_string}"
                            writer.writerow([string])
                    else:
                        if isinstance(value, list):
                            var_string = "--".join(str(x) for x in value)
                        else:
                            var_string = str(value)
                        string = f"{comment} {key}: {var_string}"
                        writer.writerow([string])

                # Write DataFrame to the file without adding an empty row
                new_df.to_csv(csvfile, index=False)

        else:
            new_df.to_csv(target_file, index=False)


def gather_buffer_files(file_location):
    # Get a list of all CSV files in the file location
    csv_files = [f for f in os.listdir(file_location) if f.endswith(".csv")]

    # Create a dictionary to store the buffer files for each original file
    buffer_files = {}
    original_file = ""
    # For each file, check if there is a corresponding buffer file
    for file in csv_files:
        if "buffer" in file.split("_"):
            # if file.endswith('_buffer.csv'):
            buffer_file = file
            original_file = file.split("_buffer_")[0] + ".csv"
            original_file_path = os.path.join(file_location, original_file)
            buffer_file_path = os.path.join(file_location, buffer_file)

            # Add the buffer file to the dictionary for the corresponding original file
            if original_file not in buffer_files:
                buffer_files[original_file] = [buffer_file_path]
            else:
                buffer_files[original_file].append(buffer_file_path)

    return original_file, buffer_files


def merge_buffer_files(file_location, check_compatibility=True):
    original_file, buffer_files = gather_buffer_files(file_location)

    # For each original file, load the original file and all buffer files into DataFrames
    for original_file, buffer_files_list in buffer_files.items():

        original_file_path = os.path.join(file_location, original_file)
        original_df = pd.read_csv(original_file_path)
        buffer_dfs = [pd.read_csv(buffer_file) for buffer_file in buffer_files_list]

        # Append the buffer file data to the original file
        for buffer_df, buffer_file in zip(buffer_dfs, buffer_files_list):
            # Check each file if it is compatible with the original file
            compatible = compare_dataframes(original_df, buffer_df)

            if compatible:
                combined_df = original_df.append(buffer_df, ignore_index=True)
            else:
                if check_compatibility:
                    print(
                        f"Warning: Buffer file {buffer_file} is inconsistent with original file. Pass check_compatibility = False to ignore this. Only use this option if you know what you're doing!"
                    )
                    return
                else:
                    combined_df = original_df.append(buffer_df, ignore_index=True)

        # Save the combined DataFrame to the original file and delete all buffer files
        combined_df.to_csv(original_file_path, index=False)
        for buffer_file in buffer_files_list:
            os.remove(buffer_file)


def test():
    df1 = pd.DataFrame(
        {
            "name": ["John", "Emily", "Jack", "Sara"],
            "age": [25, 30, 18, 22],
            "country": ["USA", "Canada", "USA", "Australia"],
        }
    )

    df2 = pd.DataFrame(
        {
            "name": ["Tom", "Kate", "Mark", "Sophie"],
            "age": [20, 35, 27, 19],
            "country": ["USA", "UK", "Canada", "France"],
        }
    )

    meta1 = {"a": [2, 3], "b": [34, 76]}
    meta2 = {"log": "True", "check": "False"}

    target_file1 = "/home/os18o068/Documents/PHD/Projects/pics_external_data/training_sets/file1.csv"
    write_to_csv(df1, target_file1, meta=meta1)

    write_to_csv(df2, target_file1, meta=meta2, conflict="add")


class FileSystem:
    def __init__(self, location, kind="basic", specs={}):
        pass


class Projects(FileSystem):
    def __init__(self, location):
        pass
