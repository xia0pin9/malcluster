# -*- coding: utf-8 -*-

import sys
cimport cython
import numpy as np
cimport numpy as np
np.import_array()
from pefile import *


cdef class ImpHash:

    cdef np.ndarray bits
    cdef int m

    def __cinit__(self, m=1):
        self.m = m
        self.bits = np.empty(256, dtype='uint8')
        for x in xrange(256):
            self.bits[x] = self.nnz(x)


    def get_imp(self, pe):
        impstrs = []
        exts = ['ocx', 'sys', 'dll']
        if not hasattr(pe, "DIRECTORY_ENTRY_IMPORT"):
            return ""
        for entry in pe.DIRECTORY_ENTRY_IMPORT:
            libname = entry.dll.lower()
            parts = libname.rsplit('.', 1)
            if len(parts) > 1 and parts[1] in exts:
                libname = parts[0]
            #print libname
            for imp in entry.imports:
                funcname = None
                if not imp.name:
                    funcname = ordlookup.ordLookup(entry.dll.lower(), imp.ordinal, make_name=True)
                    if not funcname:
                        raise Exception("Unable to look up ordinal %s:%04x" % (entry.dll, imp.ordinal))
                else:
                    funcname = imp.name

                if not funcname:
                    continue

                impstrs.append('%s.%s' % (libname.lower(),funcname.lower()))
        return impstrs


    cdef np.ndarray[np.uint8_t, ndim=1] generateHash(self, char * filename):
        cdef np.ndarray[np.uint8_t, ndim=1, mode='c'] hash = np.zeros(self.m*1024, dtype='uint8')
        cdef i, j, k, offset

        pe = PE(filename)
        impstrs_set = self.get_imp(pe)
        for k in xrange(len(impstrs_set)):
            offset = self.djb2(impstrs_set[k])
            i = offset >> 3
            j = 1 << (offset & 7)
            hash[i] |= j
            # i, j = divmod(self.djb2(ngarray[k]), 8)
            # hash[i] = hash[i] | 2**j
        return hash


    def float compareHash(self,  np.ndarray a, np.ndarray b):
        return 1 - np.sum(self.bits[a&b]) * 1.0 / np.sum(self.bits[a|b])


    cdef int djb2(self, np.ndarray[np.uint8_t, ndim=1, mode='c'] str):
        cdef int byte
        cdef int hash = 5381
        for byte in str:
            hash = (((hash << 5) + hash) + byte)  # bit array size: 1k
        return hash & (self.m*1024 - 1)


    cdef int nnz(self, unsigned char val):
        """ Count the number of 1-bits in val, val must be non-negative. """
        cdef int res = 0
        while val:
            res += 1
            val = val & (val -1)
        return res