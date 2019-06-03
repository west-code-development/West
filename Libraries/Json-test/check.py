from __future__ import print_function
import json

with open('input_file.json',"r") as file:
    idata = json.load(file)

with open('output_file.json',"r") as file:
    odata = json.load(file)

assert( idata == odata )

print("Test passed.")
