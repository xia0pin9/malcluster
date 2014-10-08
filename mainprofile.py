# -*- coding: utf-8 -*-

"""malcluster: profiling test interface"""

__version__ = "0.1.0"

import os
import sys
import shutil
import numpy as np
import mvhash
import nghash
import sdhash
import bshash
import imphash
import rlhash

import cProfile


p1 = None
p2 = None
malorder = []
fingerprints = {}
cwd = os.getcwd()
myhash = None
algorithms = [
              bshash.BsHash(81920, 7),   # BsHash works on original whole samples
              nghash.NgHash(7),          # NgHash works on original whole samples
              imphash.ImpHash(1),           # ImpHash works on original whole samples
              rlhash.RlHash(16807, 256, 1),
              mvhash.MvHash(512, 20, 0.7),  # MvHash works on python-extracted code secquences
              sdhash.SdHash()               # SdHash works on python-extracted code secquences
             ]
hash_names = ["bshash", "nghash", "imphash", "rlhash", "mvhash", "sdhash"]
imp_set = {}

def hash_gen():
    """ Wrapper for generating fingerprints from malware samples """
    shutil.rmtree(os.path.join(cwd, "hashs/"), ignore_errors=True)
    os.mkdir(os.path.join(cwd, "hashs"))
    mallist = (x for x in os.listdir(os.path.join(cwd, "samples/")))
    global imp_set
    for malware in mallist:
        malpath = os.path.join(cwd, "samples/" + malware)
        malorder.append(malware)
        fingerprints[malware] = myhash.generateRaw(malpath)
        for imp in fingerprints[malware]:
            if imp not in imp_set:
                imp_set[imp] = 1
            else:
                imp_set[imp] += 1


def main():
    global myhash 
    global imp_set

    hash_name = sys.argv[1]
    
    if hash_name in hash_names:
        myhash = algorithms[hash_names.index(hash_name)]
    else:
        print "Unknown hash name specified"
    hash_gen()
    
    for i in range(50):
        if len(imp_set) ==0:
            break
        index = max(imp_set, key=imp_set.get)
        print i, index, imp_set[index] 
        imp_set.pop(index, None)
#     aveset = []
#     for mal in fingerprints:
#         if fingerprints[mal][1] > 4000:
#             print mal, fingerprints[mal]
#         aveset.append(fingerprints[mal][4])
# #     test = fingerprints.values()
#     print "Ave: ", np.sum(aveset)/len(aveset)
#     for mal in fingerprints:
#         if not str(fingerprints[mal][2]).isdigit():
#             print mal, fingerprints[mal]
#     a = os.path.join(os.getcwd(), "samples/bettersurf-001d106b8f74e8a92d8d19b39fb4afea")
#     b = os.path.join(os.getcwd(), "samples/fosniw-fc2498ac9c1785d56441b4ac79ac2e95")
#     mvhasha = myhash.generateRaw(a)
#     mvhashb = myhash.generateRaw(b)
#     print mvhasha[:10]
#     print mvhashb[:10]
#     print myhash.compareRaw(mvhasha, mvhashb)

if __name__ == "__main__":
#     cProfile.run('main()')
    main()
