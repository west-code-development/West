import json
import sys
import os.path

def main():

    fname = str(sys.argv[1])
    key = str(sys.argv[2])

    if os.path.isfile(fname):
       with open(fname,"r") as file:
          data = json.load(file)
    else :
       print("Cannot find FILE: ",fname)
    print(data[key])

if __name__ == "__main__":
    main()

