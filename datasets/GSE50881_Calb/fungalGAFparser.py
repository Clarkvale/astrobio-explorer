# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 16:41:37 2020

@author: Benja
"""
from typing import List
import pandas as pd
import re

class ParsingError(Exception):
    pass

class GAF:
    
    def __init__(self, pathToFile: str, taxid: int):
        self.pathToFile = pathToFile
        self.taxid = taxid
        self.regex_a22 = re.compile("(C\w_\w+)")
        self.regex_sys = re.compile("(orf19\.\d+)")
        
        data = self._parse(self.pathToFile, self.taxid)
        self.df = pd.DataFrame(data)        
        
       
        
    def _tokenize(self, line: str, delim:str):
        line_dict = {}
        if delim == "\t":
            line = line.split(delim)
            line_dict["gene_name"] = line[2]
            line_dict["go_id"] = line[4]
            line_dict["description"] = line[10]
            line_dict["taxid"] = line[12]
            
        elif delim == "|":
            
            line_dict["a22_systemic_id"] = self.regex_a22.findall(line)[0]
            try:
                line_dict["systemic_id"] = self.regex_sys.findall(line)[0]
            except IndexError:
                raise ParsingError

        else:
            raise Exception("Unsupported delimiter")
        return line_dict
        
    def _parse(self, pathToFile, taxid: str):
        data_list = []
        
        with open(pathToFile, 'r') as fi:
            line = fi.readline()
            
            while line:
                if not line.startswith("!"):
                    
                    tab_tokens = self._tokenize(line, "\t")
                    
                    if tab_tokens["taxid"] != f"taxon:{taxid}":
                        line = fi.readline()
                    else:
                        try:
                             desc_tokens = self._tokenize(tab_tokens["description"], "|")
                             del tab_tokens["description"]
                             cols_out = {**tab_tokens, **desc_tokens}
                             data_list.append(cols_out)
                             line = fi.readline()
                        except ParsingError:
                            line = fi.readline()
                else:
                   line = fi.readline()
        return data_list
                    
if __name__ == "__main__":
    gaf = GAF("gene_association.cgd", 237561)
    gaf.df.to_csv("Candida_albicans_GO.csv")
             
        