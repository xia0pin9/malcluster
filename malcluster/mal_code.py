#!/usr/bin/env python

import os
import sys
import pefile
import numpy as np

def main():
    malpath = sys.argv[1]
    mallist = os.listdir(malpath)
    for mal in mallist:
        mpath = os.path.join(malpath, mal)
        pe = None
        try:
            count = 0
            output = os.path.join("samples_code", mal)
            pe = pefile.PE(mpath)
            #secnames = [s.Name.rstrip("\x00") for s in pe.sections]
            for sec in pe.sections:
                secname = sec.Name.rstrip("\x00")
                if ".text" == secname:
                    count += 1
                    malbytes = np.array(map(ord, sec.get_data()))
                    malbytes.astype('uint8').tofile(output)
                    #testbytes = None 
                    #with open(output) as f:
                    #testbytes = np.fromfile(f, dtype="uint8")
                    #print malbytes == testbytes
                    #break
                if "CODE" == secname:
                    count += 1
                    malbytes = np.array(map(ord, sec.get_data()))
                    malbytes.astype('uint8').tofile(output)
                    #testbytes = None 
                    #with open(output) as f:
                    # 	testbytes = np.fromfile(f, dtype="uint8")
                    #print malbytes == testbytes
                    #break
            if count != 1:
                print mpath, "\t", secname
        except:
            continue

if __name__ == "__main__":
    main()
