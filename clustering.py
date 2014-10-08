# -*- coding: utf-8 -*-

"""malcluster: malware clustering analysis tool"""

__version__ = "0.1.0"

import os
import sys
import time
import timeit
import shutil
import itertools
import numpy as np
import fastcluster
import multiprocessing as mp
from scipy.cluster.hierarchy import fcluster  # @UnresolvedImport
from scipy.interpolate import PiecewisePolynomial  # @UnresolvedImport
from scipy.optimize import fsolve  # @UnresolvedImport
import matplotlib.pyplot as plt
import mvhash
import nghash
import sdhash
import bshash
import imphash
import rlhash
import profile


p1 = None
p2 = None
malorder = []
fingerprints = {}
cwd = os.getcwd()
myhash = None
algorithms = [
              bshash.BsHash(81920, 10),   # BsHash works on original whole samples
              nghash.NgHash(7),          # NgHash works on original whole samples
              imphash.ImpHash(1),           # ImpHash works on original whole samples
              rlhash.RlHash(16807, 256, 1),
              mvhash.MvHash(512, 20, 0.7),  # MvHash works on python-extracted code secquences
              sdhash.SdHash()               # SdHash works on python-extracted code secquences
             ]
hash_names = ["bshash", "nghash", "imphash", "rlhash", "mvhash", "sdhash"]


def timing(f):
    def wrap(*args):
        time1 = time.time()
        ret = f(*args)
        time2 = time.time()
        print '%s function took %0.3f ms' % (f.func_name, (time2-time1)*1000.0)
        return ret
    return wrap


def hash_gen():
    """ Wrapper for generating fingerprints from malware samples """
    shutil.rmtree(os.path.join(cwd, "hashs/"), ignore_errors=True)
    os.mkdir(os.path.join(cwd, "hashs"))
    mallist = (x for x in os.listdir(os.path.join(cwd, "samples/")))
    for malware in mallist:
        malpath = os.path.join(cwd, "samples/" + malware)
        malorder.append(malware)
        fingerprints[malware] = myhash.generateHash(malpath)


def hash_comp(ab):
    return myhash.compareHash(ab[0], ab[1])


def get_dmentry():
    for item in itertools.combinations(malorder, 2):
        yield map(fingerprints.get, item)
       

def get_dmlist():
    """ Get all the pairwise distance matrix """
    number_per_round = 10000
    result = []
    getdm = get_dmentry()
    pool = mp.Pool(processes=mp.cpu_count())
    while True:
        tempresult = pool.map(hash_comp, itertools.islice(getdm, number_per_round))
        if tempresult:
            result.extend(tempresult)
        else:
            break
    return np.array(result)


def hacluster(y):
    """ Wrapper for the Hierarchical Clustering algorithm from fastcluster """
    z = fastcluster.linkage(y, method='single')
    return z

def calculate_precision(ref_lines, cdb_lines):
    tp = 0
    t_num = 0

    for cline in cdb_lines:
        clist = (cline.split(':')[1]).split()
        clist_num = len(clist)
        t_num += clist_num
        maxcnt = 0
        for rline in ref_lines:
            rlist = (rline.split(':')[1]).split()
            #rlist_num = len(rlist)
            cnt = 0
            for v in clist:
                if v in rlist:
                    cnt += 1
            if cnt > maxcnt:
                maxcnt = cnt
                if maxcnt == clist_num:
                    break
        tp += maxcnt

    return float(tp)/t_num

def calculate_recall(ref_lines, cdb_lines):
    tp = 0
    t_num = 0

    for rline in ref_lines:
        rlist = (rline.split(':')[1]).split()
        rlist_num = len(rlist)
        t_num += rlist_num
        maxcnt = 0
        for cline in cdb_lines:
            clist = (cline.split(':')[1]).split()
            #clist_num = len(clist)
            cnt = 0
            for v in rlist:
                if v in clist:
                    cnt += 1
            if cnt > maxcnt:
                maxcnt = cnt
                if maxcnt == rlist_num:
                    break
        tp += maxcnt

    return float(tp)/t_num

def evaluate(z):
    shutil.rmtree(os.path.join(os.getcwd(), "eval/"), ignore_errors=True)
    os.mkdir(os.path.join(os.getcwd(), "eval"))
    thresholds = np.linspace(0, 1, 21)
    mid = np.arange(0, 1, 0.01)
    precision_set = []
    recall_set = []
    refset = {}
    tempx = 0
    for i in xrange(len(malorder)):
        if malorder[i].split("-")[0] not in refset:
            refset[malorder[i].split("-")[0]] = [i]
        else:
            refset[malorder[i].split("-")[0]].append(i)
    with open("eval/refset.txt", "w") as f:
        for family in refset:
            f.write(family + ": " + ' '.join([str(x) for x in refset[family]]) + "\n")
    
    for i in thresholds:
        with open("eval/refset.txt") as f:
            reflines = f.readlines()
        hc = fcluster(z, i, 'distance')
        cdblines = get_clist(hc, i)
        precision = calculate_precision(reflines,  cdblines)
        recall = calculate_recall(reflines, cdblines)
        precision_set.append(precision)
        recall_set.append(recall)
    global p1
    p1 = PiecewisePolynomial(thresholds, np.array(precision_set)[:, np.newaxis])
    global p2
    p2 = PiecewisePolynomial(thresholds, np.array(recall_set)[:, np.newaxis])
    for x in mid:
        root, infodict, ier, mesg = fsolve(pdiff, x, full_output=True)
        root = float("%.3f" % root[0])
        if ier == 1 and thresholds.min() < root < thresholds.max():
            tempx = root
            break
    tempy = p2(tempx)
    if p1(tempx) > p2(tempx):
        tempy = p1(tempx)
    print "Best x:", tempx, "Best y:", tempy
    try:
        pleg1, = plt.plot(thresholds, precision_set, '-bo', label="Precision")
        pleg2, = plt.plot(thresholds, recall_set, '-rx', label="Recall")
        plt.legend([pleg1, pleg2], ["Precision", "Recall"], loc="center right")
        plt.xlabel('Distance threshold (t)')
        plt.ylabel("Precision and Recall")
        plt.show()
    except:
        raise
    finally:
        plt.close()
    return tempx, tempy


def pdiff(x):
    return p1(x) - p2(x)


def get_clist(hc, s):
    ncluster = 0
    shres = '%1.3f' % s
    with open("eval/"+shres+".txt", "w") as f:
        for i in xrange(len(malorder)):
            cid = hc[i]
            if cid == -1:
                continue
            if cid == 0:
                hc[i] = -1
                f.write("C" + str(ncluster) + ": " + str(i) + "\n")
                ncluster += 1
            else:
                f.write("C" + str(ncluster) + ": ")
                for j in range(i, len(malorder)):
                    if hc[j] == cid:
                        hc[j] = -1
                        f.write(" " + str(j))
                f.write("\n")
                ncluster += 1
    print "threshold: ", s, "Nclusters: ", ncluster
    with open("eval/"+shres+".txt") as f:
        return f.readlines()


def get_clusters(z, shresholdx):
    hc = fcluster(z, shresholdx, 'distance')
    clusters = {}
    with open("eval/cluster_results.txt", "w") as f:
        for i in xrange(hc.min(), hc.max()+1):
            clusters[i] = []
            for j in xrange(hc.size):
                if hc[j] == i:
                    mal = malorder[j]
                    clusters[i].append(mal)
        for i in clusters:
            f.write("C"+str(i)+":"+"\n")
            for mal in clusters[i]:
                f.write(mal+"\n")
            f.write("\n")
            

def main():
    global myhash 
    hash_name = sys.argv[1]
    if hash_name in hash_names:
        myhash = algorithms[hash_names.index(hash_name)]
    else:
        print "Unknown hash name provided."

    startedat = timeit.default_timer()
    hash_gen()
    hashgenat = timeit.default_timer()
    print "Finish generating fingerprints", hashgenat - startedat
    y = get_dmlist()
    print "Max distance", np.amax(y)
    getdmat = timeit.default_timer()
    print "Finish generating distance matrix", getdmat - hashgenat
    z = hacluster(y)
    hclustat = timeit.default_timer()
    shresholdx, tempy = evaluate(z)
    get_clusters(z, shresholdx)
    print "Finish clustering analyais", len(malorder), hclustat - getdmat    


if __name__ == "__main__":
#     profile.run('main()')
    main()
