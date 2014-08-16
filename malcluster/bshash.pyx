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

cdef class BsHash:

    cdef int fpSize
    cdef int shredSize
    cdef int windowSize
    # Prepare a bit number table for fast query
    cdef np.ndarray bits
    cdef float alpha


    def __cinit__(self, int fpSize = 8*1024, int shredSize = 20, int windowSize = 1):
        self.fpSize = fpSize
        self.shredSize = shredSize
        self.windowSize = windowSize
        self.bits = np.empty(256, dtype='uint8')
        for x in xrange(256):
            self.bits[x] = self.nnz(x)
    
    
    def generateHash(self, char * filename):
        """ Generate mvhash fingerprint from input file by calling the executable """
        cdef np.ndarray mvarray, mvoutput
#         inputname = os.path.join(os.getcwd(), "samples/" + filename)
#         outputname = os.path.join(os.getcwd(), "hashs/" + filename)
        outputname = filename.replace("samples", "hashs")
        commands = ["./bshash", "-s", str(self.fpSize), "-n", str(self.shredSize), "-w", str(self.windowSize), "-o", outputname, "-i", filename]
        try: 
            p = subprocess.Popen(commands, stdout=subprocess.PIPE)
            output= p.communicate()[0]     
            if len(output) != 0:
                print filename, output
        except ValueError:
            print "mvhash fingerprint generation error."
            sys.exit(1)
#         return outputname
        with open(outputname) as f:
            mvarray = np.fromfile(f, dtype="uint8")
#         if mvarray.size%256 != 0:
#             print "Invalid size of fingerprints"
#             sys.exit(1)  
#         mvoutput = np.empty(mvarray.size/256, dtype=object)    
#         for i in xrange(mvarray.size/256):
#             tempbitarray = bitarray()
#             for byte in xrange(256):
#                 tempbitarray.extend(np.binary_repr(mvarray[i*256+byte], 8)) 
#             mvoutput[i] = tempbitarray
        mvarray = mvarray[4:]
        return mvarray
      
    
    cpdef float compareHashRaw(self, char *a, char *b):
        """ Compare two bshash fingerprints through original algorithm by calling the executable """
#         fileA = os.path.join(os.getcwd(), "hashs/" + a) 
#         fileB = os.path.join(os.getcwd(), "hashs/" + b)
        commands = ["./bshash", "-i", a, "-c", b]
        try: 
            p = subprocess.Popen(commands, stdout=subprocess.PIPE)
            output= p.communicate()[0]     
            if len(output) != 0:
                return 1 - float(output)
        except:
            print "mvhash fingerprint generation error.", a, b
            sys.exit(1)
       
       
    @cython.boundscheck(False)
    cpdef float compareHash(self, a, b):
        #print np.sum(self.bits[a|b]), np.sum(self.bits[a&b]), (np.sum(self.bits[a]) + np.sum(self.bits[b]) - np.sum(self.bits[a&b]))
        return 1 - np.sum(self.bits[a&b])*1.0/ np.sum(self.bits[a|b])


    cdef int nnz(self, unsigned char val):
        cdef int res = 0
        while val:
            res += 1
            val = val & (val -1)
        return res
