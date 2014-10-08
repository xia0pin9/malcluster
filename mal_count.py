#!/usr/bin/evn python

import os
import sys

def main():
    familycount = 0
    malcount = 0
    results = {}
    mallist = os.listdir(sys.argv[1])
    for mal in mallist:
        if mal.split("-")[0] not in results:
            results[mal.split("-")[0]] = [mal]
        else:
            results[mal.split("-")[0]].append(mal)
    for family in results:
        familycount += 1
        malcount += len(results[family])
        print family.ljust(25, " "), "\t", len(results[family]) 
    print "Done, family number:", familycount, malcount

if __name__ == "__main__":
    main()
