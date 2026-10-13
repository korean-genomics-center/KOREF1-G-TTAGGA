#!/usr/bin/env python

import sys

# First arg is vcf file

rat_cut = float(sys.argv[2])

with open(sys.argv[1], "r") as f:
    for line in f:
        line = line.rstrip()
        if line.startswith("#"):
            print(line)
        else:
            fields = line.split("\t")
            tags = fields[7].split(";")

            data_format = fields[8].split(':')
            data = fields[9].split(':')

            dict_data = dict(zip(data_format, data))
            # Get supporting read info
            gt = dict_data["GT"]
            rs = dict_data["DR"]
            vs = dict_data["DV"]

            if gt == "0/0" or gt == "./.":
                continue

            rs, vs = int(rs), int(vs)
            
            rat = (vs/(rs+vs))
            if rat >= rat_cut:
                print(line)
