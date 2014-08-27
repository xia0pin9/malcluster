# -*- coding: utf-8 -*-

import sys
cimport cython
import numpy as np
cimport numpy as np
np.import_array()
import pefile


cdef class RlHash:

    cdef np.ndarray bits
    cdef unsigned int base_multiplier       # base_multiplier value
    cdef unsigned int rollingsize          # Size of n-gram featurs
    cdef unsigned int mask       # Evaluation cratiria for determine cutoffs
    cdef float fpsize           # Fingerprint size in Kbits
    cdef unsigned long modulo 

    def __cinit__(self, int base=16807, int rsize=256, float fpsize=1):
        '''
        Multiplier and modulo choosen based on experiments and reference of:
        http://en.wikipedia.org/wiki/Linear_congruential_generator
        '''
        self.base_multiplier = base
        self.rollingsize = rsize
        self.modulo = 2 ** 31 - 1
        self.fpsize = fpsize
        self.mask = self.rollingsize - 1      # Compare the lowest 11-bits, expect a average window length of 2048 bytes
        self.bits = np.empty(256, dtype='uint8')
        for x in xrange(256):
            self.bits[x] = self.nnz(x)

    
    cdef np.ndarray get_cutoffs(self, np.ndarray malbytes):
        '''
        Compute hash value of each n-gram with the following rule:
            S(k+1) = (S(k) - c_1*(b**k) + c_(k+1)) * b
        '''
        cdef unsigned long hs, base_n
        cdef unsigned int first_value, last_index, last_value, num_rolling
        
        cutoffs = []
        
        if malbytes.size == 0:
            return np.array([])
        
        if self.rollingsize < malbytes.size:
            last_index = self.rollingsize
            num_rolling = malbytes.size - self.rollingsize + 1
        else:
            last_index = malbytes.size
            num_rolling = 1

        hs = self.get_basic_rolling(malbytes[:last_index])
        # Check whether the lowest 11-bit equals to 0
        if hs & self.mask == 0:
            cutoffs.append(self.rollingsize - 1)

        first_value = malbytes[0]
        base_n = (self.base_multiplier**self.rollingsize) & self.modulo

        for i in xrange(1, num_rolling):
            last_index = i + self.rollingsize
            last_value = malbytes[last_index - 1]
            # Update current rolling hash with previous value
            hs = ((hs - first_value * base_n + last_value) * self.base_multiplier) & self.modulo
            # Check whether the lowest 11-bit equals to 0
            if hs & self.mask == 0:
                cutoffs.append(last_index)
            first_value = malbytes[i]

        # Add the rest as a single window, and mark the last byte as the last cutoff
        if last_index != malbytes.size or cutoffs == []:
            cutoffs.append(malbytes.size)

        return np.array(cutoffs, dtype='int64')


    cpdef int get_basic_rolling(self, np.ndarray str):
        '''
        Compute the first n-gram hash with the following rule:
            S(0) = c_1*(b**n) + c_2*(b**(n-1) + ... + c_k*b 
        '''
        cdef length = str.size
        cdef np.ndarray[np.int64_t, ndim=1, mode='c'] basearray  # = np.empty(length)

        basearray = np.array([(self.base_multiplier**x) & self.modulo for x in xrange(length, 0, -1)])

        return np.sum(basearray*str)  & self.modulo


    def generateHash(self, char * filename):
        cdef np.ndarray[np.uint8_t, ndim=1, mode='c'] hash = np.zeros(self.fpsize*1024, dtype='uint8')
        cdef np.ndarray[np.uint8_t, ndim=1, mode='c'] malbytes = np.array([], dtype='uint8')
        cdef np.ndarray[np.int64_t, ndim=1, mode='c'] cutoff_results = np.array([], dtype='int64')
        cdef int i, j, offset, firstpos, currentpos
        try:
            pe = pefile.PE(filename, fast_load=True)
   
            for sec in pe.sections:
                secname = sec.Name.rstrip("\x00")
                if ".text" == secname or "CODE" == secname:
                    malbytes = np.append(malbytes, np.array(map(ord, sec.get_data()[:sec.Misc_VirtualSize])).astype('uint8'))
                else:
                    continue                 
                
                cutoff_results = np.append(cutoff_results, self.get_cutoffs(malbytes))
                firstpos = 0
                for currentpos in cutoff_results:
                    offset = self.djb2(malbytes[firstpos:currentpos])
                    i = offset >> 3
                    j = 1 << (offset & 7)
                    hash[i] |= j
                    firstpos = currentpos
        except:
            raise
        
        return hash


    cpdef float compareHash(self,  np.ndarray a, np.ndarray b):
#         if np.sum(self.bits[a|b]) == 0:
#             return 1
#         else:
#             return 1 - np.sum(self.bits[a&b]) * 1.0 / np.sum(self.bits[a|b])
        cdef int bitwise_and = np.sum(self.bits[a&b])
        #print np.sum(self.bits[a|b]), np.sum(self.bits[a&b]), (np.sum(self.bits[a]) + np.sum(self.bits[b]) - np.sum(self.bits[a&b]))
        return float("%.6f" % (bitwise_and*1.0/ (np.sum(self.bits[a]) + np.sum(self.bits[b]) - bitwise_and)))



    cdef int djb2(self, np.ndarray str):
        cdef int byte
        cdef int hash = 5381
        for byte in str:
            hash = (((hash << 5) + hash) + byte)  # bit array size: 1k
        return hash & (int)(self.fpsize*1024 - 1)


    cdef int djb2a(self, np.ndarray str):
        cdef int b
        cdef int h = 5381
        for b in str:
            h = (((h << 5) ^ h) ^ b)  # bit array size: 32*1024*8
        return h & (int)(self.fpsize*1024-1)


    cdef int nnz(self, unsigned char val):
        """ Count the number of 1-bits in val, val must be non-negative. """
        cdef int res = 0
        while val:
            res += 1
            val = val & (val -1)
        return res
