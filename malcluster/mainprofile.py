# -*- coding: utf-8 -*-

"""malcluster: malware clustering analysis tool"""

__version__ = "0.1.0"

import os
import shutil
import numpy as np
import mvhash
import nghash
import sdhash
import bshash
import imphash
import rlhash

import cProfile
import pstats


p1 = None
p2 = None
malorder = []
fingerprints = {}
cwd = os.getcwd()
myhash = None


def hash_gen():
    """ Wrapper for generating fingerprints from malware samples """
    shutil.rmtree(os.path.join(cwd, "profile_hashs/"), ignore_errors=True)
    os.mkdir(os.path.join(cwd, "profile_hashs"))
    mallist = (x for x in os.listdir(os.path.join(cwd, "profile_samples/")))
    for malware in mallist:
        malpath = os.path.join(cwd, "profile_samples/" + malware)
        malorder.append(malware)
        fingerprints[malware] = myhash.generateHash(malpath)


def main():
    global myhash 
#     myhash = mvhash.MvHash(512, 20, 0.7)
    myhash = nghash.NgHash(16)
#     myhash = sdhash.SdHash()
#     myhash = bshash.BsHash(81920, 16)
#     myhash = imphash.ImpHash(1)
#     myhash = rlhash.RlHash(16807, 2048, 2)
#     hash_gen()
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
    a = os.path.join(os.getcwd(), "samples/bettersurf-001d106b8f74e8a92d8d19b39fb4afea")
    b = os.path.join(os.getcwd(), "samples/fosniw-fc2498ac9c1785d56441b4ac79ac2e95")
    mvhasha = myhash.generateHash(a)
    mvhashb = myhash.generateHash(b)
    print myhash.compareHash(mvhasha, mvhashb)

if __name__ == "__main__":
    cProfile.run('main()')
#     main()
