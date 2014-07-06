# -*- coding: utf-8 -*-

"""malcluster: malware clustering analysis tool"""

__version__ = "0.1.0"

import os
import time
import timeit
import shutil
import itertools
import numpy as np
import fastcluster
import quality
import multiprocessing as mp
from scipy.cluster.hierarchy import fcluster
from scipy.interpolate import PiecewisePolynomial
from scipy.optimize import fsolve 
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist
from matplotlib.pyplot import figure, show
import mvhash
import nghash
import sdhash
import bshash

def timing(self, f):
    def wrap(*args):
        time1 = time.time()
        ret = f(*args)
        time2 = time.time()
        print '%s function took %0.3f ms' % (f.func_name, (time2-time1)*1000.0)
        return ret
    return wrap

p1 = None
p2 = None
malorder = []
fingerprints = {}
cwd = os.getcwd()
myhash = None
    
#     def __init__(self):
#         pass
#         if hashname == 'mvhash':
#             hash = MvHash()
#         elif hashname == 'nghash':
#             hash = NgHash()

""" Wrapper for generating fingerprints from malware samples """
def hashGen():
    shutil.rmtree(os.path.join(cwd, "hashs/"), ignore_errors=True)
    os.mkdir(os.path.join(cwd, "hashs"))
    mallist = (x for x in os.listdir(os.path.join(cwd,"samples/")))
    for malware in mallist:
        malpath = os.path.join(cwd, "samples/" + malware)
        malorder.append(malware)
        fingerprints[malware] = myhash.generateHash(malpath)

def hashComp(ab):
#         a = fingerprints[ab[0]]
#         b = fingerprints[ab[1]]
    return myhash.compareHash(ab[0], ab[1])
    
def getDmentry():
    for item in itertools.combinations(malorder, 2):
        yield map(fingerprints.get, item)
       
""" Get all the pairwise distance matrix """
def getDmlist():
    numberPerRound = 10000
    result = []
    getdm = getDmentry()
    pool = mp.Pool(processes=mp.cpu_count())
    while True:
        tempresult = pool.map(hashComp, itertools.islice(getdm, numberPerRound))
        if tempresult:
            result.extend(tempresult)
        else:
            break
    return np.array(result)

# def getDmlist():
#     dmlist = [map(fingerprints.get, item) for item in itertools.combinations(malorder, 2)]
#     pool = mp.Pool(processes=mp.cpu_count())
#     result = pool.map(hashComp, dmlist)
#     pool.close()
#     pool.join()
#     print "DM size: ", len(result)
#     return np.array(result)

""" Wrapper for the Hierarchical Clustering algorithm from fastcluster """
def hacluster(y):
    Z = fastcluster.linkage(y, method='complete')
    return Z

def evaluate(Z):
    shutil.rmtree(os.path.join(os.getcwd(), "eval/"), ignore_errors=True)
    os.mkdir(os.path.join(os.getcwd(), "eval"))
    threshold = np.linspace(0, 1, 21)
    mid = np.arange(0, 1, 0.01)
    precisionSet = []
    recallSet = []
    refset = {}
    reflines = None
    tempx = 0
    for i in xrange(len(malorder)):
        if malorder[i].split("-")[0] not in refset:
            refset[malorder[i].split("-")[0]] = [i]
        else:
            refset[malorder[i].split("-")[0]].append(i)
    with open("eval/refset.txt", "w") as f:
        for family in refset:
            f.write(family + ": " + ' '.join([str(x) for x in refset[family]]) + "\n")
    
    for i in threshold:
        cdblines = None
        with open("eval/refset.txt") as f:
            reflines = f.readlines()
        hc = fcluster(Z, i, 'distance')
        cdblines = getClist(hc, i)
        precision = quality.precision(reflines,  cdblines)
        recall = quality.recall(reflines, cdblines)
#         precision, recall = clusterStats(np.array(hc))
        precisionSet.append(precision)
        recallSet.append(recall)
#     figure = plt.figure()
    global p1
    p1 = PiecewisePolynomial(threshold, np.array(precisionSet)[:, np.newaxis])
    global p2
    p2 = PiecewisePolynomial(threshold, np.array(recallSet)[:, np.newaxis])
    for x in mid:
        root, infodict, ier, mesg = fsolve(pdiff, x, full_output=True)
        root = float("%.3f" % root[0])
        if ier == 1 and threshold.min() < root < threshold.max():
            tempx = root
            break
    tempy = p2(tempx)
    if p1(tempx) > p2(tempx):
        tempy = p1(tempx)
    print "Best x:", tempx, "Best y:", tempy
    pleg1, = plt.plot(threshold, precisionSet, '-bo', label="Precision")
    pleg2, = plt.plot(threshold, recallSet, '-rx', label="Recall")
    plt.legend([pleg1, pleg2], ["Precision", "Recall"], loc="center right")
    plt.xlabel('Distance Threshold (t)')
    plt.ylabel("Precision and Recall")
    plt.show()
    plt.savefig('myfig.pdf')
    return tempx , tempy
#     print precisionSet, recallSet

def pdiff(x):
    return p1(x) - p2(x)

def getClist(hc, s):
    ncluster = 0
    shres = '%1.3f' % s
    with open("eval/"+shres+".txt", "w") as f:
        for i in xrange(len(malorder)):
            cid = hc[i]
            if (cid == -1):
                continue
            if (cid == 0):
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
    print "Threshold: ", s, "Nclusters: ", ncluster
    with open("eval/"+shres+".txt") as f:
        return f.readlines()
            
def getClusters(Z, shresholdx):
    hc = fcluster(Z, shresholdx, 'distance')
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

# def clusterStats(hc):
#     global malorder
#     tcluster = {}
#     ccluster = {}
#     precision = recall = 0
#     for mal in malorder:
#         if mal.split("-")[0] not in tcluster:
#             tcluster[mal.split("-")[0]] = [mal.split("-")[1]]
#         else:
#             tcluster[mal.split("-")[0]].append(mal.split("-")[1])
#     for i in xrange(hc.min(), hc.max()):
#         ccluster[i] = {}
#         for j in xrange(hc.size):
#             if hc[j] == i:
#                 mal = malorder[j]
#                 if mal.split("-")[0] not in ccluster[i]:
#                     ccluster[i][mal.split("-")[0]] = [mal.split("-")[1]]
#                 else:
#                     ccluster[i][mal.split("-")[0]].append(mal.split("-")[1])
#     for i in ccluster:
#         tempNum = 0
#         tempName = None
#         for mals in ccluster[i]:
#             if len(ccluster[i][mals]) > tempNum:
#                 tempNum = len(ccluster[i][mals])
#                 tempName = mals
#         if len(tcluster[tempName]) > tempNum:
#             precision += tempNum
#         else:
#             print "Got error when calculating precision value for", tempName
#     for malName in tcluster:
#         tempNum = 0
#         for j in ccluster:
#             if malName in ccluster[j]:
#                 if len(ccluster[j][malName]) > tempNum:
#                     tempNum = len(ccluster[j][malName])
#         if tempNum < len(tcluster[malName]):
#             recall += tempNum
#         else:
#             print "Got error when calculating recall value for", malName
#     return (precision*1.0/len(malorder), recall*1.0/len(malorder))
            

def main():
    global myhash 
#     myhash = mvhash.MvHash(512, 20, 0.7)
#     myhash = nghash.NgHash(20, 80)
#     myhash = sdhash.SdHash()
    myhash = bshash.BsHash(81920, 20, 1)
    startedat = timeit.default_timer()
    hashGen()
    hashgenat = timeit.default_timer()
    print "Finish generating fingerprints", hashgenat - startedat
    y = getDmlist()
    print max(y)
    getdmat = timeit.default_timer()
    print "Finish generating distance matrix", getdmat - hashgenat
    Z = hacluster(y)
    hclustat = timeit.default_timer()
    shresholdx, tempy = evaluate(Z)
    #print "Round", i, j, tempy
    getClusters(Z, shresholdx)  
    print "Finish clustering analyais", len(malorder), hclustat - getdmat    
    
#     currentij = ''
#     besty = 0
#     for i in range(1000, 4000):         
#         for j in range(5, 55):
#             try:
#                 myhash = mvhash.MvHash(i, j, 0.7)
# #     myhash = nghash.NgHash(16, 8)
# #     myhash = sdhash.SdHash()
#     #startedat = timeit.default_timer()
#                 hashGen()
#     #hashgenat = timeit.default_timer()
#     #print "Finish generating fingerprints", hashgenat - startedat
#                 y = getDmlist()
#     #print max(y)
#     #getdmat = timeit.default_timer()
#     #print "Finish generating distance matrix", getdmat - hashgenat
#                 Z = hacluster(y)
#     #hclustat = timeit.default_timer()
#                 shresholdx, tempy = evaluate(Z)
#                 print "Round", i, j, tempy
#     #getClusters(Z, shresholdx)
#                 if tempy > besty:
#                     besty = tempy
#                     currentij = str(i) + "-" + str(j)
#             except:
#                 continue	
#     print "Done", currentij, besty

#     print "Test1"
#     b = os.path.join(os.getcwd(), "samples/nsis-2.45-setup.exe")
#     a = os.path.join(os.getcwd(), "samples/nsis-2.46-setup.exe")
#     mvhasha = myhash.generateHash(a)
#     mvhashb = myhash.generateHash(b)
#     print myhash.compareHashH(mvhasha, mvhashb)
#     nghash = NgHash(16, 32)
#     print "Test2"
#     rawa = nghash.generateRaw(a)
#     rawb = nghash.generateRaw(b)
#     print nghash.compareRaw(rawa, rawb)
#     print "Test3"
#     hasha = nghash.generateHash(a)
#     hashb = nghash.generateHash(b)
#     print nghash.compareHashH(hasha, hashb)
#     print "Test4"
#     hasha = nghash.generateHasha(a)
#     hashb = nghash.generateHasha(b)
#     print nghash.compareHashJ(hasha, hashb)


if __name__ == "__main__":
    main()
