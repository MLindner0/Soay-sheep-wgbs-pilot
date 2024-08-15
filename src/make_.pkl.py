"""
Short script to turn the .csv file created in R (see pipeline-data_from_filenames.R) into a .pkl object for easier processing within the snakemake pipeline
"""

# import libraries
import csv
import pickle

# read . csv file, format data, and save .pkl file
new_data_dict = {}
with open("/shared/slate_group1/Shared/methylated_soay/soay_wgbs_pilot_mar2023/src/ReadGroups.help.new.csv", 'r') as data_file:
    data = csv.DictReader(data_file, delimiter=",")
    for row in data:
        item = new_data_dict.get(row["Sample"], dict())
        item[row["Key"]] = row["Value"]
        new_data_dict[row["Sample"]] = item

print(new_data_dict)
dict_ReadGroups = new_data_dict

# save
with open('src/dict_ReadGroups.pkl', 'wb') as fp:
    pickle.dump(dict_ReadGroups, fp)
    print('dictionary saved successfully to file')
