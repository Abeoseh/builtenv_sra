"""Fixes the taxonomy file from qiime2"""

## note: if you run this 2x or more you'll get an error since the files that are made are placed in the same directory as the unprocessed files.
import pandas as pd
import glob
import re

files = glob.glob("./csv_files/combine/*.txt")
print(f"found files: {files}")

for file in files:
    df = pd.read_csv(file, sep="\t", skiprows = 1)

    # split taxnonmy into multiple columns
    split_taxonomy = ["Domain", "Phylum", "Class","Order" ,"Family", "Genus","Species" ]
    df[split_taxonomy] = df["taxonomy"].str.split(";", expand=True)
    
    # remove taxnonmy column
    df.drop("taxonomy", axis=1, inplace=True)    
    # print(df.head())

    # remove d__, p__, etc.
    for column in split_taxonomy:
        df[column] = (df[column].str.split("__", expand=True)[1])   

    # remove chloroplast and other non-bacteria
    split_taxonomy =  [taxonomy for taxonomy in split_taxonomy if taxonomy != "Genus"]
    df = df[(df.Domain.isin(["Bacteria"]))]
    df = df[ ~df["Genus"].isin(["Mitochondria", "Chloroplast"]) ]
    df = df.dropna(subset=["Genus"])

    df.drop(split_taxonomy, axis = 1, inplace=True)

    # remove OTU ID column
    df.pop("#OTU ID")

    # rename columns so all the files don't have the same file names 
    filename = file.split("/")[3].split(".")[0] 
    genus_column = df.pop("Genus")
    print(filename)
    
    # rename columns 
    rename = {}
    with open(f"./{filename}/manifest.tsv", "r") as file:
        for line in file:

            line = line.strip("\n")
            line = line.split("\t")
            line[1] = re.sub(r".+\/", "", line[1])
            line[1] = re.sub(r"\..+", "", line[1])
            line[1] = re.sub(r"\_.+", "", line[1])

            rename[line[0]] = line[1]

    rename.pop("sample-id")
    df = df.rename(columns=rename)


    # move genus to the front
    df.insert(0,"Genus", genus_column)

    # sum common genera
    df = df.groupby("Genus", as_index=False).sum()
    df = df.reset_index()

    df.pop("index")
    
    df.to_csv(f'./csv_files/combine/single_{filename}.csv', index=False)

    print(f"done with file: {filename}")

print("done with all.")

