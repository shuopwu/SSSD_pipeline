import pandas as pd
import numpy as np
import os
import sys

def blocks(files, size=65536):
    while True:
        b = files.read(size)
        if not b: break
        yield b

def count_rows(files):
    d = {'file':[], 'how_many_rows':[]}
    for file in files:
        with open(file, "r") as f:
            d['file'].append(os.path.basename(file))

            count = sum(bl.count("\n") for bl in blocks(f)) - 1
            d['how_many_rows'].append(count) 
            
    return pd.DataFrame(d)

if __name__ == "__main__":
    count_rows(file_list)