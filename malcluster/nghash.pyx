# -*- coding: utf-8 -*-

import sys
import os
cimport cython
import subprocess
import numpy as np
cimport numpy as np
from bitarray import bitarray
np.import_array()

cdef class NgHash:

    # n-gram size, eg 16-gram
    cdef int shredSize
    cdef int windowSize
    # Prepare a bit number table for fast query
    cdef np.ndarray bits

    @cython.boundscheck(False)
    def __cinit__(self, int shredSize = 20, int windowSize = 1):
        self.shredSize = shredSize
        self.windowSize = windowSize
        self.bits = np.empty(256, dtype='uint8')
        cdef int x
        for x in xrange(256):
            self.bits[x] = self.nnz(x)

    cdef np.ndarray[np.uint8_t, ndim=1] getBytes(self, char * filename):
#         inputname = os.path.join(os.getcwd(), "samples/" + filename)
        cdef np.ndarray[np.uint8_t, ndim=1] filebytes
        with open(filename) as f:
            filebytes = np.fromfile(f, dtype="uint8")
        return filebytes
        
    def generateRaw(self, char * filename):
        cdef np.ndarray filebytes = self.getBytes(filename)
        cdef int numNg = filebytes.size - self.n + 1
        cdef np.ndarray ngarray = np.empty(numNg, dtype='S'+str(2*self.n))
        for i in xrange(numNg):
            ngarray[i] = ''.join(format(x, '02x') for x in filebytes[i:i+self.n])
        return np.unique(ngarray)
    
    @cython.boundscheck(False)    
    cpdef float compareRaw(self, np.ndarray a, np.ndarray b):
        return 1 - np.intersect1d(a, b, assume_unique=True).size *1.0/np.union1d(a, b).size    
    
    
    @cython.boundscheck(False)
    def generateHashb(self, char * filename):  
        """ Generate feature hashing fingerprint through dhb2()"""         
        cdef np.ndarray[np.uint8_t, ndim=1, mode='c'] hash = np.zeros(self.m*1024, dtype='uint8')
        cdef np.ndarray[np.uint8_t, ndim=1, mode='c'] filebytes = self.getBytes(filename)
        cdef int i, j, k, offset, numNg = filebytes.size - self.n + 1
        cdef np.ndarray[np.uint8_t, ndim=2, mode='c'] ngarray = np.empty((numNg, self.n), dtype='uint8')
#         cdef np.ndarray[np.uint8_t, ndim=1, mode='c'] ngram = np.empty(self.n, dtype='uint8')
        for i in xrange(numNg):
            ngarray[i] = filebytes[i:i+self.n]
        #print filebytes.size, ngarray.size, ngarray[0], ngarray[1][0]
        for k in xrange(numNg):
#             offset = self.djb2(ngarray[k])
#             i = offset >> 3
#             j = 1 << (offset & 7)
#             hash[i] |= j
            i, j = divmod(self.djb2(ngarray[k]), 8)
            hash[i] = hash[i] | 2**j  
        return hash
    
    
    @cython.boundscheck(False)
    def generateHasha(self, char * filename):
        """ Generate feature hashing fingerprint through djb2a()"""
        cdef np.ndarray[np.uint8_t, ndim=1, mode='c'] filebytes = self.getBytes(filename)   
        cdef int i, k, j, l, numBlock, numNg = filebytes.size - self.n + 1
        cdef np.ndarray[np.uint8_t, ndim=1, mode='c'] hasha # = np.zeros(self.m*1024, dtype='uint8')
        cdef np.ndarray[np.uint8_t, ndim=2, mode='c'] ngarraya = np.empty((numNg, self.n), dtype='uint8')
        cdef np.ndarray mvoutput
        if numNg % 1000 == 0:
            numBlock = numNg / 1000
        else:
            numBlock = numNg / 1000 + 1
            
        hasha = np.zeros(numBlock*2048, dtype='uint8')
        
#         print "Test"
        for i in xrange(numNg):
            ngarraya[i] = filebytes[i:i+self.n]
            
        for k in xrange(numNg):
#             hasha[i + l*2048] |= j            
            i, j = divmod(self.djb2a(ngarraya[k]), 8)
            l = k / 1000
            hasha[i + l*2048] = hasha[i + l*2048] | 2**j
            
        mvoutput = np.empty(hasha.size/256, dtype=object)    
        for i in xrange(hasha.size/256):
            tempbitarray = bitarray()
            for byte in xrange(256):
                tempbitarray.extend(np.binary_repr(hasha[i*256+byte], 8)) 
            mvoutput[i] = tempbitarray
        return mvoutput
    
    
    def generateHash(self, char * filename):
        """ Generate mvhash fingerprint from input file by calling the executable """
        cdef np.ndarray mvarray, mvoutput
#         inputname = os.path.join(os.getcwd(), "samples/" + filename)
#         outputname = os.path.join(os.getcwd(), "hashs/" + filename)
        outputname = filename.replace("samples", "hashs")
        commands = ["./nghash", "-n", str(self.shredSize), "-w", str(self.windowSize), "-o", outputname, "-i", filename]
        try: 
            p = subprocess.Popen(commands, stdout=subprocess.PIPE)
            output= p.communicate()[0]     
            if len(output) != 0:
                print filename, output
        except ValueError:
            print "nghash fingerprint generation error."
            sys.exit(1)
        except:
            raise
        return outputname
#         with open(outputname) as f:
#             ngarray = np.fromfile(f, dtype="uint8")
#         if mvarray.size%256 != 0:
#             print "Invalid size of fingerprints"
#             sys.exit(1)  
#         mvoutput = np.empty(mvarray.size/256, dtype=object)    
#         for i in xrange(mvarray.size/256):
#             tempbitarray = bitarray()
#             for byte in xrange(256):
#                 tempbitarray.extend(np.binary_repr(mvarray[i*256+byte], 8)) 
#             mvoutput[i] = tempbitarray
        return outputname


    cpdef float compareHash(self, char *a, char *b):
        """ Compare two nghash fingerprints through original algorithm by calling the executable """
        commands = ["./nghash", "-i", a, "-c", b]
        try: 
            p = subprocess.Popen(commands, stdout=subprocess.PIPE)
            output= p.communicate()[0]     
            if len(output) != 0:
                return float(output)
        except:
            print "mvhash fingerprint generation error.", a, b
            sys.exit(1)
      
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.nonecheck(False)       
    cpdef float compareHasha(self, np.ndarray a,
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

    @cython.boundscheck(False)
    cpdef float compareHashb(self, a, b):
        #print np.sum(self.bits[a|b]), np.sum(self.bits[a&b]), (np.sum(self.bits[a]) + np.sum(self.bits[b]) - np.sum(self.bits[a&b]))
        return 1 - np.sum(self.bits[a&b])*1.0/ np.sum(self.bits[a|b])
    
    @cython.boundscheck(False)
    cpdef float compareHashH(self, a, b):
        return np.sum(self.bits[a^b])*1.0/np.sum(self.m*1024*8)
        
    cdef int nnz(self, unsigned char val):
        cdef int res = 0
        while val:
            res += 1
            val = val & (val -1)
        return res
        
    cdef int djb2(self, np.ndarray[np.uint8_t, ndim=1, mode='c'] str):
        cdef int byte
        cdef int hash = 5381
        for byte in str:
            hash = (((hash << 5) + hash) + byte)  # bit array size: 32*1024*8
        return hash & (self.m*1024*8-1)

    cdef int djb2a(self, np.ndarray[np.uint8_t, ndim=1, mode='c'] str):
        cdef int b
        cdef int h = 5381
        for b in str:
            h = (((h << 5) ^ h) ^ b) % (2048) # bit array size: 32*1024*8
        return h
