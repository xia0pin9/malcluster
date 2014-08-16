# -*- coding: utf-8 -*-

#from __future__ import division
import subprocess
import sys
import os
import cython
from bitarray import bitarray
import numpy as np
cimport numpy as np
np.import_array()

cdef class MvHash:

    cdef int numPerBF
    cdef int neighborSize
    # Prepare a bit number table for fast query
    cdef np.ndarray bits
    cdef float alpha

    def __cinit__(self, int numPerBF = 1536, int neighborSize = 20, float alpha = 0.5):
        self.numPerBF = numPerBF
        self.neighborSize = neighborSize
        self.alpha = alpha
        self.bits = np.empty(256, dtype='uint8')
        for x in xrange(256):
            self.bits[x] = self.nnz(x)
    
    """ Generate mvhash fingerprint from input file by calling the executable """
    def generateHash(self, char * filename):
        cdef np.ndarray mvarray, mvoutput
#         inputname = os.path.join(os.getcwd(), "samples/" + filename)
#         outputname = os.path.join(os.getcwd(), "hashs/" + filename)
        outputname = filename.replace("samples", "hashs")
        commands = ["./mvhash", "-e", str(self.numPerBF), "-w", str(self.neighborSize), "-t", "8", "-o", outputname, "-i", filename]
        try: 
            p = subprocess.Popen(commands, stdout=subprocess.PIPE)
            output= p.communicate()[0]     
            if len(output) != 0:
                print filename, output
        except ValueError:
            print "mvhash fingerprint generation error."
            sys.exit(1)
        with open(outputname) as f:
            mvarray = np.fromfile(f, dtype="uint8")
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
      
    """ Compare two mvhash fingerprints through original algorithm by calling the executable """
    cpdef int compare1(self, char *a, char *b):
#         fileA = os.path.join(os.getcwd(), "hashs/" + a) 
#         fileB = os.path.join(os.getcwd(), "hashs/" + b)
        commands = ["./mvhash", "-e", str(self.numPerBF), "-i", a, "-m", b]
        try: 
            p = subprocess.Popen(commands, stdout=subprocess.PIPE)
            output= p.communicate()[0]     
            if len(output) != 0:
                if int(output) == -1:
                    return 100
                else:
                    return 100 - int(output)
        except:
            print "mvhash fingerprint generation error.", a, b
            sys.exit(1)
       
    """ Compare two mvhash fingerprints through hamming distance, normalized to 1 """   
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.nonecheck(False)       
    cpdef float compareHash(self, np.ndarray a,
                                   np.ndarray b):
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
            #return float("%.3f" % (hammingTotal/nblocksB * ratio + 1 - ratio)) # *self.alpha + (1-ratio)*(1-self.alpha)))  #
            return float("%.3f" % ((hammingTotal + (nblocksA-nblocksB))/nblocksA))
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
            #return float("%.3f" % (hammingTotal/nblocksA * ratio + 1 - ratio)) # *self.alpha + (1-ratio)*(1-self.alpha)))  #
            return float("%.3f" % ((hammingTotal + (nblocksB-nblocksA))/nblocksB) )
            
    """ Compare Two mvhash fingerprints through bit-wise jaccard similarity, normalized to 1 """
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.nonecheck(False)
    cpdef float compareHashJ(self, np.ndarray a,
                                   np.ndarray b):

        cdef float ratio, jaccardTotal = 0
        cdef int numArrayA = a.size, numArrayB = b.size
        cdef int j, i, k, index
        cdef np.ndarray[np.float_t, ndim=2] resultHamming
    
        if (numArrayA > numArrayB):
            ratio = numArrayB*1.0/numArrayA
            resultHamming = np.empty((numArrayB, numArrayA), dtype=np.float)
            for i in xrange(numArrayB):
                for j in xrange(numArrayA):
                    resultHamming[i][j] = (b[i]&a[j]).count(True)*1.0/(b[i]|a[j]).count(True)
            # Recursively find k smallest distance from a nblocksB by nblocksA array
            for k in xrange(numArrayB):
                index = np.argmin(resultHamming)
                i, j = divmod(index, numArrayA)
#                 i, j = divmod(index, numArrayA-k)
                jaccardTotal += resultHamming[i][j]
#                 resultHamming = np.delete(np.delete(resultHamming, i, 0), j, 1)
                resultHamming[i, :] = 1
                resultHamming[:, j] = 1
#             if np.round(1 - jaccardTotal/numArrayB, 4) + ratio > 1:
#                 return 1
#             else:
            return float( "%.3f" % (jaccardTotal/numArrayA)) # * ratio) #+ 1 - ratio))
        else:
            ratio = numArrayA*1.0/numArrayB
            resultHamming = np.empty((numArrayA, numArrayB), dtype=np.float)
            for i in xrange(numArrayA):
                for j in xrange(numArrayB):
                    resultHamming[i][j] = (a[i]&b[j]).count(True)*1.0/(a[i]|b[j]).count(True)
            # Recursively find k smallest distance from a nblocksA by nblocksB array
            for k in xrange(numArrayA):
                index = np.argmin(resultHamming)
                i, j = divmod(index, numArrayB)
#                 i, j = divmod(index, numArrayB-k)
                jaccardTotal += resultHamming[i][j]
#                 resultHamming = np.delete(np.delete(resultHamming, i, 0), j, 1)
                resultHamming[i, :] = 1
                resultHamming[:, j] = 1
#             if np.round(1 - jaccardTotal/numArrayA, 4) + ratio > 1:
#                 return 1
#             else:
            return float("%.3f" % (jaccardTotal/numArrayB)) #* ratio) # + 1 - ratio))

    """ Count the number of 1-bits in val, val must be non-negative. """
    cdef int nnz(self, unsigned char val):
        cdef int res = 0
        while val:
            res += 1
            val = val & (val -1)
        return res
