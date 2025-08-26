import pandas as pd
import os, time
import subprocess
from zlib_ng import gzip_ng


def run_fastp_report(list_files, path_save):
    os.makedirs(os.path.dirname(path_save), exist_ok=True)
    with gzip_ng.open(path_save, "wb") as fwb:
        for file in list_files:
            with gzip_ng.open(file, "rb") as frb:
                for line in frb:
                    fwb.write(line)


if __name__ == "__main__":
    import argparse
    
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument("--files", nargs='+')
    argument_parser.add_argument("--out")
    
    args = argument_parser.parse_args()
    
    run_fastp_report(args.files, args.out)