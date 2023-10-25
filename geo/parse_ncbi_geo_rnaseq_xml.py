#!/usr/bin/env python3

import xml.etree.ElementTree as ET
import pandas as pd
import argparse
import os

parse = parser = argparse.ArgumentParser(
    description="Extract sample information from NCBI GEO database RNA SEQ family.xml file ",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)

parse.add_argument(
    "-x", "--xml", type=str, required=True, help="path to family.xml file "
)

parse.add_argument(
    "-d", "--outdir", type=str, required=True, help="path to result csv file  "
)

parse.add_argument(
    "-p", "--prefix", type=str, required=True, help="path to result csv file  "
)

args = parse.parse_args()

xml = os.path.abspath(args.xml)
id = args.prefix
out_dir = os.path.abspath(args.outdir)

tree = ET.parse(xml)
root = tree.getroot()

info = dict()
for sample in root.iter("{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Sample"):
    info[sample.attrib["iid"]] = dict()
    for child in sample:
        if child.tag.split("}")[1] == "Title":
            # print(child.text)
            title = child.text
            info[sample.attrib["iid"]]["title"] = title
        for j in child:
            if j.tag.split("}")[1] == "Characteristics":
                # print(j,j.tag,j.text.strip(),j.items()[0][1],sample.attrib["iid"])
                # if j.items()[0][1]=="cohort":
                info[sample.attrib["iid"]][j.items()[0][1]] = j.text.strip()

result = pd.DataFrame.from_dict(info)
result = result.T
result.to_csv(f"{out_dir}/{id}_info.csv")
