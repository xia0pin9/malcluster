# -*- coding: utf-8 -*-

import subprocess
import sys
import cython
import base64
from bitarray import bitarray
import numpy as np
cimport numpy as np
np.import_array()

cdef class SdHash:

    def generateHash(self, char * filename):
        outputname = filename.replace("samples", "hashs")
        commands = ["sdhash", filename, "-o", outputname]
        try:
            p = subprocess.Popen(commands, stdout = subprocess.PIPE)
            output = p.communicate()[0]
            if len(output) != 0:
                print filename, output
            #return outputname + ".sdbf"
        except:
            raise
            #print "sdhash fingerprint generation error."
            #sys.exit(1)
        with open(outputname+".sdbf") as f:
            content = f.read().rstrip().split(":")[-1]
            mvarray = np.array(map(ord, base64.b64decode(content)))
        if mvarray.size%256 != 0:
            print "Invalid size of fingerprints"
            sys.exit(1)  
        mvoutput = np.empty(mvarray.size/256, dtype=object)    
        for i in xrange(mvarray.size/256):
            tempbitarray = bitarray()
            for byte in xrange(256):
                tempbitarray.extend(np.binary_repr(mvarray[i*256+byte], 8)) 
            mvoutput[i] = tempbitarray
        return mvoutput


    cpdef float compareHashJ(self, char *a, char *b):
        commands = ["sdhash", "-c", a, b, "-t", "0"]
        try:
            p = subprocess.Popen(commands, stdout = subprocess.PIPE)
            output = p.communicate()[0]
            score = 0
            if len(output) != 0:
                score = output.rstrip().split("|")[2]
                score = int(score[0])*100 + int(score[1])*10 + int(score[2])
            return float(score*1.0/100)
        except:
            raise
            #print "sdhash fingerprint comparison error.", a, b
            #sys.exit(1)
            

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.nonecheck(False)       
    cpdef float compareHash(self, np.ndarray a,
                                   np.ndarray b):
        """ Compare two mvhash fingerprints through hamming distance, normalized to 1 """
        cdef float ratio, hammingTotal = 0
        cdef int nblocksA = a.size, nblocksB = b.size
        cdef int i, j, k
        cdef np.ndarray[np.float_t, ndim=2] resultHamming   
    
        if (nblocksA > nblocksB):
            ratio = nblocksB*1.0/nblocksA
            resultHamming = np.empty((nblocksB, nblocksA), dtype=np.float)
            for i in xrange(nblocksB):
                for j in xrange(nblocksA):
                    resultHamming[i][j] = (b[i]^a[j]).count(True)/ 2048.0
            # Recursively find k smallest distance from a nblocksB by nblocksA array
            for k in xrange(nblocksB):
#                 print nblocksA-k, np.shape(resultHamming)[1]
                index = np.argmin(resultHamming)
                i, j = divmod(index, nblocksA)
#                 i, j = divmod(index, nblocksA-k)
                hammingTotal += resultHamming[i][j]
#                 resultHamming = np.delete(np.delete(resultHamming, i, 0), j, 1)
                resultHamming[i, :] = 1
                resultHamming[:, j] = 1
#             if np.round(hammingTotal/nblocksB, 3) * self.alpha + ratio * (1) > 1:
#                 return 1
#             else:
            return float("%.3f" % (hammingTotal/nblocksB * ratio + 1 - ratio)) # *self.alpha + (1-ratio)*(1-self.alpha)))  #
        else:
            ratio = nblocksA*1.0/nblocksB
            resultHamming = np.empty((nblocksA, nblocksB), dtype=np.float)
            for i in xrange(nblocksA):
                for j in xrange(nblocksB):
                    resultHamming[i][j] = (a[i]^b[j]).count(True)/ 2048.0
            # Recursively find k smallest distance from a nblocksA by nblocksB array
            for k in xrange(nblocksA):
                index = np.argmin(resultHamming)
                i, j = divmod(index, nblocksB)
#                 i, j = divmod(index, np.shape(resultHamming)[1])
                hammingTotal += resultHamming[i][j]
#                 resultHamming = np.delete(np.delete(resultHamming, i, 0), j, 1)
                resultHamming[i, :] = 1
                resultHamming[:, j] = 1
#             if np.round(hammingTotal/nblocksA, 3)*ratio + 1 - ratio > 1:
#                 return 1
#             else:
            return float("%.3f" % (hammingTotal/nblocksA  * ratio + 1 - ratio)) # *self.alpha + (1-ratio)*(1-self.alpha)))  #