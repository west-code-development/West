#!/usr/bin/python3

import json
import yaml
import datetime

def dump_message(fileName,docText,meta={}) :
    data = {}
    try : 
       data = yaml.load(docText,Loader=yaml.SafeLoader) 
    except :
       print("Cannot interpret docText :",docText)
    # append
    with open(fileName, 'a+') as file:
       file.write("---\n")
       file.write("time: "+datetime.datetime.utcnow().replace(microsecond=0).isoformat()+"\n")
       file.write("meta: "+json.dumps(meta)+"\n")
       file.write("data: "+json.dumps(data)+"\n")

def append_log(*args, **kwargs):
    #
    fileName = args[0]
    jsonText = args[1]
    if "meta" in kwargs.keys() : 
       dump_message(fileName,jsonText,meta=kwargs["meta"])
    else : 
       dump_message(fileName,jsonText)

def clear_file(*args, **kwargs):
    #
    fileName = args[0]
    with open(fileName, 'w') as file:
       pass

def test() :
    #
    fileName = "west.log.yaml"
    jsonText = '{ "a" : "b", "c" : [1,2], "d" : { "e": 2, "f" : [3.5677,3.909090900], "l" : true } }'
    clear_file(fileName)
    append_log(fileName,jsonText)
    append_log(fileName,jsonText,meta={"type" : "setup" })
    append_log(fileName,jsonText,meta={"type" : "setup", "ciao" : True })
    append_log(fileName,jsonText,meta={"ciao" : False })

if __name__ == "__main__":
    # execute only if run as a script
    test()

