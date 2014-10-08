# -*- coding: utf-8 -*-

import sys
cimport cython
import numpy as np
cimport numpy as np
from bitarray import bitarray
np.import_array()


cdef class NgHash:

    # n-gram size, eg 16-gram
    cdef int ng_size
    cdef int fp_size
    cdef float fpbitsrate 
    # Prepare a bit number table for fast query
    cdef np.ndarray bits


    @cython.boundscheck(False)
    def __cinit__(self, int ngsize = 16):
        self.ng_size = ngsize
        self.fp_size = 2048                 # Default single block length in bits
        self.fpbitsrate = 0.7
        self.bits = np.empty(256, dtype='uint8')
        cdef int x
        for x in xrange(256):
            self.bits[x] = self.nnz(x)


    cdef np.ndarray[np.uint8_t, ndim=1] get_bytes(self, char * filename):
        cdef np.ndarray[np.uint8_t, ndim=1, mode='c'] filebytes
        with open(filename) as f:
            filebytes = np.fromfile(f, dtype="uint8")
        return filebytes


    def generateRaw(self, char * filename):
        cdef np.ndarray filebytes = self.get_bytes(filename)
        cdef int num_ng = filebytes.size - self.ng_size + 1
        cdef np.ndarray ngarray = np.empty(num_ng, dtype='S'+str(2*self.ng_size))
        cdef int i 
 
        for i in xrange(num_ng):
            ngarray[i] = ''.join(format(x, '02x') for x in filebytes[i:i+self.ng_size])
        return np.unique(ngarray)
    
   
    cpdef float compareRaw(self, np.ndarray a, np.ndarray b):
        return 1 - np.intersect1d(a, b, assume_unique=True).size *1.0/np.union1d(a, b).size    
    
    
    def generateHash(self, char * filename):  
        """ Generate feature hashing fingerprint through djb2()"""         
#         cdef np.ndarray[np.uint8_t, ndim=1, mode='c'] hash
        #cdef np.ndarray[np.uint8_t, ndim=1, mode='c'] malbytes = np.array([], dtype='uint8')
        cdef np.ndarray malbytes = self.get_bytes(filename)
        cdef np.ndarray ng_hash = np.array([])
        cdef np.ndarray ngoutput
        cdef np.ndarray[np.uint8_t, ndim=1, mode='c'] temp_block
#         cdef np.ndarray[np.uint8_t, ndim=1, mode='c'] filebytes = self.get_bytes(filename)
        cdef int i, j, k, offset, byte, num_ng, firstpos = 0
        cdef np.ndarray[np.uint8_t, ndim=2, mode='c'] ngarray
        cdef int block_index, block_count, fp_count
        
        try:
            num_ng = malbytes.size - self.ng_size + 1  
            ngarray = np.empty((num_ng, self.ng_size), dtype='uint8')
          
            for i in xrange(num_ng):
                ngarray[i] = malbytes[i:i+self.ng_size]

            fp_count = (int)(self.fp_size * self.fpbitsrate*1.0 + 0.5)
            block_count = num_ng / fp_count
            if num_ng % fp_count != 0:
                block_count += 1

            for block_index in xrange(block_count):
                temp_block = np.zeros(self.fp_size/8, dtype='uint8')
                for k in xrange(block_index*fp_count, (block_index + 1)*fp_count):
                    if k == num_ng:
                        break
                    else:
                        offset = self.djb2(ngarray[k])
                        i = offset >> 3
                        j = 1 << (offset & 7)
                        temp_block[i] |= j
#                 firstpos = (block_index + 1)*fp_count
                ng_hash = np.concatenate([ng_hash, temp_block])

            if ng_hash.size%256 != 0:
                print "Invalid size of fingerprints"
                sys.exit(1)  
            
            ngoutput = np.empty(ng_hash.size/256, dtype=object)    
            for i in xrange(ng_hash.size/256):
                tempbitarray = bitarray()
                for byte in xrange(256):
                    tempbitarray.extend(np.binary_repr(ng_hash[i*256+byte], 8)) 
                ngoutput[i] = tempbitarray
        except:
            raise

        return ngoutput

     
    cpdef float compareHash(self, np.ndarray a, np.ndarray b):
        """ Compare two mvhash fingerprints through hamming distance, normalized to 1 """
        cdef float ratio, hammingTotal = 0
        cdef int nblocksA = a.size, nblocksB = b.size
        cdef int i, j, k
        cdef np.ndarray[np.float_t, ndim=2] resultHamming
        #print "Test compare hash size: ", a.size, b.size
        if (nblocksA > nblocksB):
            ratio = nblocksB*1.0/nblocksA
            resultHamming = np.empty((nblocksB, nblocksA), dtype=np.float)
            for i in xrange(nblocksB):
                for j in xrange(nblocksA):
                    resultHamming[i][j] = (b[i]^a[j]).count(True)/ (self.fp_size*1.0)
            # Recursively find k smallest distance from a nblocksB by nblocksA array
            for k in xrange(nblocksB):
                index = np.argmin(resultHamming)
                #print "Test compare hash index value: ", index
                i, j = divmod(index, nblocksA)
                hammingTotal += resultHamming[i][j]
                resultHamming[i, :] = 1
                resultHamming[:, j] = 1
            return float("%.8f" % ((hammingTotal + (nblocksA-nblocksB))/nblocksA))
        else:
            ratio = nblocksA*1.0/nblocksB
            resultHamming = np.empty((nblocksA, nblocksB), dtype=np.float)
            for i in xrange(nblocksA):
                for j in xrange(nblocksB):
                    resultHamming[i][j] = (a[i]^b[j]).count(True)/ (self.fp_size*1.0)
            # Recursively find k smallest distance from a nblocksA by nblocksB array
            for k in xrange(nblocksA):
                index = np.argmin(resultHamming)
                #print "Test compare hash index value: ", index
                i, j = divmod(index, nblocksB)
                hammingTotal += resultHamming[i][j]
                resultHamming[i, :] = 1
                resultHamming[:, j] = 1
            return float("%.6f" % ((hammingTotal + (nblocksB-nblocksA))/nblocksB) )    


    cdef int nnz(self, unsigned char val):
        cdef int res = 0
        while val:
            res += 1
            val = val & (val -1)
        return res


    cdef int djb2(self, np.ndarray[np.uint8_t, ndim=1, mode='c'] str):
        cdef unsigned int byte
        cdef int hash = 5381
        for byte in str:
            hash = (((hash << 5) + hash) + byte)  # bit array size: 32*1024*8
        return hash & (self.fp_size-1)


    cdef int djb2a(self, np.ndarray[np.uint8_t, ndim=1, mode='c'] str):
        cdef unsigned int b
        cdef int h = 5381
        for b in str:
            h = (((h << 5) ^ h) ^ b) % (2048) # bit array size: 32*1024*8
        return h
