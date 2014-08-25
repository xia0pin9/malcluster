# -*- coding: utf-8 -*-

import sys
cimport cython
import numpy as np
cimport numpy as np
np.import_array()
from pefile import *


cdef class RlHash:

    cdef np.ndarray bits
    cdef unsigned int base       # prim base value
    cdef unsigned int n          # Size of n-gram featurs
    cdef unsigned int mask       # Evaluation cratiria for determine cutoffs
    cdef unsigned long module 

    def __cinit__(self, b):
        '''
        Multiplier and moduler choosen based on experiments and reference of:
        http://en.wikipedia.org/wiki/Linear_congruential_generator
        '''
        self.base = b
        self.module = 2 ** 31 - 1
        self.n = 2048
        self.mask = self.n - 1      # Compare the lowest 11-bits, expect a average window length of 2048 bytes
        self.bits = np.empty(256, dtype='uint8')
        for x in xrange(256):
            self.bits[x] = self.nnz(x)


    cdef np.ndarray[np.uint8_t, ndim=1] getBytes(self, char * filename):
        cdef np.ndarray[np.uint8_t, ndim=1] filebytes
        with open(filename) as f:
            filebytes = np.fromfile(f, dtype="uint8")
        return filebytes

    
    cdef np.ndarray get_cutoffs(self, np.ndarray malbytes):
        '''
        Compute hash value of each n-gram with the following rule:
            S(k+1) = (S(k) - c_1*(b**k) + c_(k+1)) * b
        '''
        cdef np.ndarray[np.uint8_t, ndim=1, mode='c'] ngarray = np.empty(self.n, dtype='uint8')
        cdef unsigned long hs,
        cdef unsigned int first, last, num_rolling
        
        cutoffs = []
        
        if malbytes.size == 0:
            return np.array([])
        
        if self.n < malbytes.size:
            last = self.n
            num_rolling = malbytes.size - self.n + 1
        else:
            last = malbytes.size
            num_rolling = 1

        hs = self.get_basic_rolling(malbytes[:last])
        # Check whether the lowest 11-bit equals to 0
        if hs & self.mask == 0:
            cutoffs.append(self.n - 1)

        first = malbytes[0]

        for i in xrange(1, num_rolling):
            ngarray = malbytes[i: i + self.n]
            last = ngarray[self.n - 1]
            # Update current rolling hash with previous value
            hs = ((hs - first * (self.base**self.n) + last) * self.base) & self.module
            # Check whether the lowest 11-bit equals to 0
            if hs & self.mask == 0:
                cutoffs.append(i + self.n - 1)
            first = ngarray[0]

        # Add the rest as a single window, and mark the last byte as the last cutoff
        if last != malbytes[malbytes.size - 1] or cutoffs == []:
            cutoffs.append(malbytes.size - 1)

        return np.array(cutoffs, dtype='int64')


    cpdef int get_basic_rolling(self, np.ndarray str):
        '''
        Compute the first n-gram hash with the following rule:
            S(0) = c_1*(b**n) + c_2*(b**(n-1) + ... + c_k*b 
        '''
        cdef length = str.size
        cdef np.ndarray[np.int64_t, ndim=1, mode='c'] basearray  # = np.empty(length)

        basearray = np.array([(self.base**x) & self.module for x in xrange(length, 0, -1)])

        return np.sum(basearray*str)  & self.module


    def generateHash(self, char * filename):
#         cdef np.ndarray[np.uint8_t, ndim=1, mode='c'] hash = np.zeros(self.m*1024, dtype='uint8')
        cdef np.ndarray[np.uint8_t, ndim=1, mode='c'] malbytes = np.array([], dtype='uint8')
        cdef np.ndarray[np.int64_t, ndim=1, mode='c'] cutoff_results
        cdef int numNg
        cdef i, j, k, offset
        try:
#             pe = PE(filename)
#   
#             for sec in pe.sections:
#                 secname = sec.Name.rstrip("\x00")
#                 if ".text" == secname:
#                     malbytes = np.append(malbytes, np.array(map(ord, sec.get_data())).astype('uint8'))
#                 if "CODE" == secname:
#                     malbytes = np.append(malbytes, np.array(map(ord, sec.get_data())).astype('uint8')) 
                                    
            malbytes = self.getBytes(filename)
            cutoff_results = self.get_cutoffs(malbytes)
            
#             cut_length = []
#             firstpos = 0
#             for i in cutoff_results:
# #                 feature = malbytes[firstpos: i]
#                 cut_length.append(i - firstpos)
#                 firstpos = i + 1
#             print "Test of sum: ", malbytes.size, sum(cut_length), cutoff_results[cutoff_results.size-1], len(cut_length)
            
        except:
            raise
        
        return cutoff_results.size


    cpdef float compareHash(self,  np.ndarray a, np.ndarray b):
        return 1
#         if np.sum(self.bits[a|b]) == 0:
#             return 1
#         else:
#             return 1 - np.sum(self.bits[a&b]) * 1.0 / np.sum(self.bits[a|b])


    cdef int djb2(self, np.ndarray str):
        cdef int byte
        cdef int hash = 5381
        for byte in str:
            hash = (((hash << 5) + hash) + byte)  # bit array size: 1k
        return hash & (self.m*1024 - 1)


    cdef int djb2a(self, np.ndarray str):
        cdef int b
        cdef int h = 5381
        for b in str:
            h = (((h << 5) ^ h) ^ b)  # bit array size: 32*1024*8
        return h & (self.m*1024-1)


    cdef int nnz(self, unsigned char val):
        """ Count the number of 1-bits in val, val must be non-negative. """
        cdef int res = 0
        while val:
            res += 1
            val = val & (val -1)
        return res
