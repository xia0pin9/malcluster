# -*- coding: utf-8 -*-

"""malcluster: malware clustering analysis tool"""

__version__ = "0.1.0"

import os
import shutil
import cPickle
import itertools
import numpy as np
import multiprocessing as mp
import matplotlib.pyplot as plt
import mvhash
import nghash
import sdhash
import bshash
import imphash
import rlhash


malorder = []
fingerprints = {}
cwd = os.getcwd()  # @UndefinedVariable
myhash = None  # @UndefinedVariable
hash_names = ["bshash", "nghash", "imphash", "rlhash", "mvhash", "sdhash"]


def hash_gen():
    """ Wrapper for generating fingerprints from malware samples
    :rtype : fingerprints
    """
    shutil.rmtree(os.path.join(cwd, "hashs/"), ignore_errors=True)  # @UndefinedVariable
    os.mkdir(os.path.join(cwd, "hashs"))
    mallist = (x for x in os.listdir(os.path.join(cwd, "samples/")))
    for malware in mallist:
        malpath = os.path.join(cwd, "samples/" + malware)
        fingerprints[malware] = myhash.generateHash(malpath)


def hash_comp(ab):
    return myhash.compareHash(ab[0], ab[1])


def get_dmentry_same():
    """
    Get dmentry for same family samples, malorder is a single family set
    :return: fingerprints pair from the same family
    """
    for item in itertools.combinations(malorder, 2):  # @UndefinedVariable
        yield map(fingerprints.get, item)


def get_dmentry_different(mal_families):
    """
    Get dmentry for different family samples
    :return: fingerprints pair from different families
    """
    for (i, j) in itertools.combinations(mal_families.keys(), 2):
        for temp_i in mal_families[i]:
            for temp_j in mal_families[j]:
                yield (fingerprints.get(temp_i), fingerprints.get(temp_j))


def get_dmlist(mal_families):
    """ Get all the pairwise distance matrix
    :rtype : distance matrics
    """
    number_per_round = 10000
    result = []
    if not mal_families:
        getdm = get_dmentry_same()
    else:
        getdm = get_dmentry_different(mal_families)
    pool = mp.Pool(processes=mp.cpu_count())
    while True:  # @UndefinedVariable
        tempresult = pool.map(hash_comp, itertools.islice(getdm, number_per_round))
        if tempresult:
            result.extend(tempresult)
        else:
            break
    return result


def generate_data(hashindex):
    global malorder
    global myhash
    algorithms = [
                bshash.BsHash(81920, 7, 1),   # BsHash works on original whole samples
                nghash.NgHash(7, 1),          # NgHash works on original whole samples
                imphash.ImpHash(1),           # ImpHash works on original whole samples
                rlhash.RlHash(16807, 256, 1),
                mvhash.MvHash(512, 20, 0.7),  # MvHash works on python-extracted code secquences
                sdhash.SdHash()               # SdHash works on python-extracted code secquences
                ]
    
    if hashindex == -1:
        hash_algos = algorithms
    else:
        hash_algos = algorithms[hashindex]
       
    # if len(os.listdir("samples_whole")) == 1146 and hash_algo:
    #     for file_name in os.listdir("samples"):
    #         shutil.move("samples/" + file_name, "samples_code/" + file_name)
    #         shutil.move("samples_whole/" + file_name, "samples/" + file_name)
        
    for hash_type in hash_algos:
        hash_name = hash_names[hash_algos.index(hash_type)]
        
#         if j == 2:
#             if len(os.listdir("samples_code")) == 1146:
#                 for file_name in os.listdir("samples"):
#                     shutil.move("samples/" + file_name, "samples_whole/" + file_name)
#                     shutil.move("samples_code/" + file_name, "samples/" + file_name)
        
        print "Generating fingerprint lists for %s (%d)." % (hash_name, len(os.listdir("samples")))  # @UndefinedVariable
        myhash = hash_type
        hash_gen()
        
        for mal in fingerprints:
            mal_family = mal.split("-")[0]
            if mal_family not in mal_families:
                mal_families[mal_family] = [mal]
            else:
                mal_families[mal_family].append(mal)
    
        # Calculate the distances within each single family
        same_family_dm = []
        for family in mal_families:
            malorder = mal_families[family]
            print "Calculating pairwise distance for family %s (%d)." % (family, len(mal_families[family]))
            same_family_dm.extend(get_dmlist(None))
        with open(hash_name + ".same", 'w+b') as f:
            cPickle.dump(same_family_dm, f)
            
        # Calculate the pairwise distances between every two families
        print "Calculating pairwise distance for samples from different families."
        diff_family_dm = get_dmlist(mal_families)
        with open(hash_name + ".diff", 'w+b') as f:
            cPickle.dump(diff_family_dm, f)


def display_original():
    for hash_name in hash_names:
        with open(hash_name + ".same", 'r+b') as f:
            same_family_dm = cPickle.load(f)
        same_family_uniqw, same_family_inverse = np.unique(same_family_dm, return_inverse=True)
        plt.figure(0)
        plt.plot(same_family_uniqw, np.bincount(same_family_inverse), label=hash_name)
        plt.legend(loc='upper right')
        plt.title("Same Family Distances Distribution")
        plt.xlabel("Distance")
        plt.ylabel("Number Count")
            
        with open(hash_name + ".diff", 'r+b') as f:
            diff_family_dm = cPickle.load(f)
        diff_family_uniqw, diff_family_inverse = np.unique(diff_family_dm, return_inverse=True)
        plt.figure(1)
        plt.plot(diff_family_uniqw, np.bincount(diff_family_inverse), label=hash_name)
        plt.legend(loc='upper left')
        plt.title("Diff Family Distances Distribution")
        plt.xlabel("Distance")
        plt.ylabel("Number Count")
    plt.show()


def display_cdf():
    for hash_name in hash_names:
        with open(hash_name + ".same", 'r+b') as f:
            same_family_dm = cPickle.load(f)
        print "Calculate same family dmcount matrics for hash", hash_name 
        dmcount_total = len(same_family_dm)
        same_family_uniqw, same_family_inverse = np.unique(same_family_dm, return_inverse=True)
        same_family_dmcount = dict(zip(same_family_uniqw, np.bincount(same_family_inverse)*1.0/dmcount_total))
      
        plt.figure(0)
        same_family_x = np.sort(np.array(same_family_dmcount.keys()))
        same_family_y = np.zeros(same_family_x.size)
        same_family_y[0] = same_family_dmcount[same_family_x[0]]
        for i in xrange(1, same_family_x.size):
            same_family_y[i] = same_family_y[i-1] + same_family_dmcount[same_family_x[i]]
        plt.plot(same_family_x, same_family_y, label=hash_name)
        plt.legend(loc='lower right')
        plt.title("Same Family Evaluation (recall)")
        plt.xlabel("Distance")
        plt.ylabel("Cumulative Probability")
 
        with open(hash_name + ".diff", 'r+b') as f:
            diff_family_dm = cPickle.load(f)
        print "Calculate diff family dmcount matrics for hash", hash_name
        dmcount_total = len(diff_family_dm)
        diff_family_uniqw, diff_family_inverse = np.unique(diff_family_dm, return_inverse=True)
        diff_family_dmcount = dict(zip(diff_family_uniqw, np.bincount(diff_family_inverse)*1.0/dmcount_total))

        plt.figure(1)
        diff_family_x = np.sort(np.array(diff_family_dmcount.keys()))
        diff_family_y = np.zeros(diff_family_x.size)
        diff_family_y[0] = diff_family_dmcount[diff_family_x[0]]
        for i in xrange(1, diff_family_x.size):
            diff_family_y[i] = diff_family_y[i-1] + diff_family_dmcount[diff_family_x[i]]
        plt.plot(diff_family_x, diff_family_y, label=hash_name)
        plt.legend(loc='upper left')
        plt.title("Different Families Evaluation (precision)")
        plt.xlabel("Distance")
        plt.ylabel("Cumulative Probability")

    plt.show()
    

def main():
    generate_data()
    display_original()
    display_cdf()

    print "Finish fuzzy hashing evaluation analyais"


if __name__ == "__main__":
    main()
