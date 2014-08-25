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


    cdef get_imp(self, pe):
        impstrs = []
#         exts = ['ocx', 'sys', 'dll', 'drv']

        for entry in getattr(pe, "DIRECTORY_ENTRY_IMPORT", []):
            libname = entry.dll.lower()
            if 'invalid' in libname:
                continue
            for imp in entry.imports:
                if not imp.name:
                    continue
                else:
                    funcname = imp.name.lower()

                impstrs.append('%s.%s' % (libname, funcname))
                
#         exps = getattr(pe, "DIRECTORY_ENTRY_EXPORT", [])    
#         if exps != []:
#             for func in exps.symbols:  
#                 if not func.name:
#                     continue
#                 else:
#                     impstrs.append(func.name.lower())  
#                  
#         for bound_imp_desc in getattr(pe, "DIRECTORY_ENTRY_BOUND_IMPORT", []):
#             libname = bound_imp_desc.name
#             for bound_imp_ref in bound_imp_desc.entries:
#                 if bound_imp_ref.name:
#                     impstrs.append('%s.%s' % (libname.lower(), bound_imp_ref.name.lower()))
#                     
#         for module in getattr(pe, "DIRECTORY_ENTRY_DELAY_IMPORT", []):
#             libname = module.dll.lower()
#             for symbol in module.imports:
#                 if symbol.name:
#                     impstrs.append('%s.%s' % (libname, symbol.name.lower()))
# 
#         resources = getattr(pe, "DIRECTORY_ENTRY_RESOURCE", [])
#         if resources != []:
#             for entry in resources.entries:
#                 if entry.name != None:
#                     impstrs.append(str(entry.name).lower())
#             directory = getattr(entry, 'directory', [])
#             if directory != []:
#                 for id in directory.entries:
#                     impstrs.append(str(id.name).lower())

        return impstrs


    def generateHash(self, char * filename):
        cdef np.ndarray[np.uint8_t, ndim=1, mode='c'] hash = np.zeros(self.m*1024, dtype='uint8')
        cdef int i, j, k, offset
        try:
            pe = PE(filename)
            impstrs_set = self.get_imp(pe)
            for k in xrange(len(impstrs_set)):
                impstr = np.array(list(bytearray(impstrs_set[k])))
                offset = self.djb2(impstr)
                i = offset >> 3
                j = 1 << (offset & 7)
                hash[i] |= j
                # i, j = divmod(self.djb2(ngarray[k]), 8)
                # hash[i] = hash[i] | 2**j
        except:
            raise
        return hash


    cpdef float compareHash(self,  np.ndarray a, np.ndarray b):
        if np.sum(self.bits[a|b]) == 0:
            return 1
        else:
            return 1 - np.sum(self.bits[a&b]) * 1.0 / np.sum(self.bits[a|b])


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
        return h & (self.m*1024 - 1)


    cdef int nnz(self, unsigned char val):
        """ Count the number of 1-bits in val, val must be non-negative. """
        cdef int res = 0
        while val:
            res += 1
            val = val & (val -1)
        return res
