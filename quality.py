#!/usr/bin/python
#
# quality.py
#  : measure clustering precision / recall
#
#  - Jiyong Jang (2012)
#

import sys
import os

def precision(ref_lines, cdb_lines):
    tp = 0
    t_num = 0

    for cline in cdb_lines:
        clist = (cline.split(':')[1]).split()
        clist_num = len(clist)
        t_num += clist_num
        maxcnt = 0
        for rline in ref_lines:
            rlist = (rline.split(':')[1]).split()
            rlist_num = len(rlist)
            cnt = 0
            for v in clist:
                if v in rlist:
                    cnt += 1
            if cnt > maxcnt:
                maxcnt = cnt
                if maxcnt == clist_num:
                    break
        tp += maxcnt

    return float(tp)/t_num

def recall(ref_lines, cdb_lines):
    tp = 0
    t_num = 0

    for rline in ref_lines:
        rlist = (rline.split(':')[1]).split()
        rlist_num = len(rlist)
        t_num += rlist_num
        maxcnt = 0
        for cline in cdb_lines:
            clist = (cline.split(':')[1]).split()
            clist_num = len(clist)
            cnt = 0
            for v in rlist:
                if v in clist:
                    cnt += 1
            if cnt > maxcnt:
                maxcnt = cnt
                if maxcnt == rlist_num:
                    break
        tp += maxcnt

    return float(tp)/t_num

