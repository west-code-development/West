#!/usr/bin/python3

#
# Copyright (C) 2015-2023 M. Govoni
# This file is distributed under the terms of the
# GNU General Public License. See the file `License'
# in the root directory of the present distribution,
# or http://www.gnu.org/copyleft/gpl.txt .
#
# This file is part of WEST.
#

import json
from xml.etree import ElementTree as ET


def jsonString2data(jsonString):
    try:
        data = json.loads(jsonString)
    except:
        print("Cannot convert jsonString to data: ", jsonString)
    return data


def function3D_to_base64(*args, **kwargs):
    #
    fileName = args[0]
    #
    data = {}
    #
    root = ET.parse(fileName)
    grid_function = root.find("grid_function")
    data["grid_function"] = grid_function.text.replace("\n", "")
    assert grid_function.attrib["type"] in ["double", "complex"]
    data["dtype"] = grid_function.attrib["type"]
    grid = root.find("grid")
    data["grid"] = [int(grid.attrib["nx"]), int(grid.attrib["ny"]), int(grid.attrib["nz"])]
    domain = root.find("domain")
    data["domain"] = {
        "a": [float(f) for f in domain.attrib["a"].split()],
        "b": [float(f) for f in domain.attrib["b"].split()],
        "c": [float(f) for f in domain.attrib["c"].split()],
    }
    #
    return data


def base64_to_function3D(*args, **kwargs):
    #
    fileName = args[0]
    #
    # root
    attrib = {}
    attrib["name"] = kwargs["name"]
    attrib["xmlns:xsi"] = "http://www.w3.org/2001/XMLSchema-instance"
    attrib["xsi:schemaLocation"] = "http://www.quantum-simulation.org/ns/fpmd/fpmd-1.0 function3d.xsd"
    root = ET.Element("{http://www.quantum-simulation.org/ns/fpmd/fpmd-1.0}function3d", attrib=attrib)
    # domain
    data = jsonString2data(kwargs["domain"])
    attrib = {}
    for l in ["a", "b", "c"]:
        attrib[l] = f"{data[l][0]} {data[l][1]} {data[l][2]}"
    ET.SubElement(root, "domain", attrib=attrib)
    # grid
    data = jsonString2data(kwargs["grid"])
    attrib = {}
    attrib["nx"] = f"{data[0]}"
    attrib["ny"] = f"{data[1]}"
    attrib["nz"] = f"{data[2]}"
    ET.SubElement(root, "grid", attrib=attrib)
    # grid_function
    attrib = {}
    assert kwargs["dtype"] in ["double", "complex"]
    attrib["type"] = kwargs["dtype"]
    attrib["nx"] = f"{data[0]}"
    attrib["ny"] = f"{data[1]}"
    attrib["nz"] = f"{data[2]}"
    attrib["encoding"] = "base64"
    ET.SubElement(root, "grid_function", attrib=attrib).text = kwargs["grid_function"]
    # write
    ET.ElementTree(root).write(fileName, encoding="UTF-8", xml_declaration=True)
    #
    return 0


def test():
    #
    base64_to_function3D(
        "vext.xml",
        name="delta_v",
        domain='{"a":[1,0,0],"b":[0,1,0],"c":[0,0,1]}',
        grid="[2,3,3]",
        grid_function="encoded\nfunction\ngoes\nhere\n",
        dtype="double",
    )
    print(function3D_to_base64("vext.xml"))


if __name__ == "__main__":
    # execute only if run as a script
    test()
