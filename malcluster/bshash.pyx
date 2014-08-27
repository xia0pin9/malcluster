# -*- coding: utf-8 -*-

#from __future__ import division
import subprocess
import sys
import os
import pefile 
import cython
from bitarray import bitarray
import numpy as np
cimport numpy as np
np.import_array()

cdef class BsHash:

    cdef int fp_size
    cdef int ng_size
    # Prepare a bit number table for fast query
    cdef np.ndarray bits
    cdef float alpha


    def __cinit__(self, int fpsize = 80*1024, int ngsize = 16):
        self.fp_size = fpsize           # Fingerprint size in bytes, 80KB
        self.ng_size = ngsize       
        self.bits = np.empty(256, dtype='uint8')
        for x in xrange(256):
            self.bits[x] = self.nnz(x)


    def generateHash(self, char * filename):
        """ Generate mvhash fingerprint from input file by calling the executable """
        cdef np.ndarray[np.uint8_t, ndim=1, mode='c'] malbytes = np.array([], dtype='uint8')
        cdef np.ndarray bshash = np.zeros(self.fp_size, dtype='uint8')
        cdef int i, j, k, offset, num_ng 

        try: 
            pe = pefile.PE(filename, fast_load=True)

            for sec in pe.sections:
                secname = sec.Name.rstrip("\x00")
                if ".text" == secname or "CODE" == secname:
                    malbytes = np.append(malbytes, np.array(map(ord, sec.get_data()[:sec.Misc_VirtualSize])).astype('uint8'))
                else:
                    continue

            num_ng = malbytes.size - self.ng_size + 1

            for k in xrange(num_ng):
                offset = self.djb2(malbytes[k:k+self.ng_size])  
                i = offset >> 3
                j = 1 << (offset & 7)
                bshash[i] |= j
        except:
            raise
        return bshash


    @cython.boundscheck(False)
    cpdef float compareHash(self, np.ndarray a, np.ndarray b):
#         cdef int bitwise_union = np.sum(self.bits[a|b])
#         
#         if bitwise_union == 0:
#             return 1
#         else:
#             return 1 - np.sum(self.bits[a&b]) * 1.0 / bitwise_union
        cdef int bitwise_and = np.sum(self.bits[a&b])
  
        return float("%.6f" % (1 - bitwise_and*1.0/ (np.sum(self.bits[a]) + np.sum(self.bits[b]) - bitwise_and)))


    cdef int djb2(self, np.ndarray[np.uint8_t, ndim=1, mode='c'] str):
        cdef unsigned int byte
        cdef int hash = 5381
        
        for byte in xrange(str.size):
            hash = (((hash << 5) + hash) + str[byte]) 
        return hash & (self.fp_size*8-1)


    cdef int nnz(self, unsigned char val):
        cdef int res = 0
        while val:
            res += 1
            val = val & (val -1)
        return res
