#!/usr/bin/python

import sys, getopt, glob, re, sqlite3
import pandas as pd



def main(argv):
    #initialzing
    inputdir = ''
    try:
       opts, args = getopt.getopt(argv,"hi:",["idir="])
    except getopt.GetoptError:
       print('meta_build.py -i <inputdir>')
       sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print("meta_build.py -i <inputdir>")
            sys.exit()
        elif opt in ("-i", "--idir"):
          inputdir = arg
    #getting files in dir
    filenames = [file for file in glob.iglob(fr"{inputdir}\*\*_meta.txt")]

    #parsing files iteravely and storing line dicts in a big list

    api_blacklist = ["Escherichia coli", "Streptococcus mutans", "Candida albicans"]
    dict_list = [parse(file) for file in filenames]
    api_available = ["no" if x["org"] in api_blacklist else "yes" for x in dict_list]

    for x in range(len(dict_list)):
        dict_list[x]["GSE_Analysis_Compatible"] = api_available[x] 


    df = pd.DataFrame(dict_list)

    makeDB(df)

    
def parse(filename):
    re_dict = {"\[1] \"ORGANISM: (\w*\s\w*)": 0, 
    "\[1] \"MICROGRAVITY TYPE:(\s\w*)": 0, 
    "type:\s*\n\s*Expression profiling by (.*)": re.M,
    "geo_accession:\s*(.*)": re.M}
    line_dict = {"org":None, "grav": None, "type": None, "acc":None}
    with open(filename, "r") as f:
        file = f.read()
        for key in zip(dict.keys(re_dict), dict.keys(line_dict)):
            reg, entry = key
            regex = re.compile(reg, flags = re_dict[reg])
            line_dict[entry] = regex.findall(file)[0]
    

    return line_dict


def makeDB(pandasDF):
    conn = sqlite3.connect("AstroMeta.db")
    c = conn.cursor()
    c.execute("DROP TABLE IF EXISTS STUDIES")
    c.execute("CREATE TABLE STUDIES (GEO_Acc text, Organism text, Grav_Type text, Exp_Type text)")
    conn.commit()
    pandasDF.to_sql("STUDIES", conn, if_exists = "replace", index = False)
    conn.close()


   
   

if __name__ == "__main__":
   main(sys.argv[1:])

    