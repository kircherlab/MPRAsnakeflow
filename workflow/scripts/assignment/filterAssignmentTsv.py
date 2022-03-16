#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

:Author: Martin Kircher
:Contact: martin.kircher@bihealth.de
:Date: *17.09.2013
"""

import sys, os
from optparse import OptionParser
from collections import defaultdict
import gzip

## Assumes a barcode sorted file with three columns
## 1: barcode 
## 2: assigned insert
## 3: quality metrics of assignment (i.e. alignment)

parser = OptionParser("%prog [options]")
parser.add_option("-m","--minimum", dest="minimum", help="Minimum support for an assignment to get reported (def 3)",default=3,type="int")
parser.add_option("-f","--fraction", dest="fraction", help="Fraction of reads required to support assignment (def 0.75)",default=0.75,type="float")
parser.add_option("-o","--other", dest="other", help="Also report barcodes for unknown/other insert sequence (def off)",default=False,action="store_true")
parser.add_option("-a","--ambiguous", dest="ambiguous", help="Also report barcodes for ambiguous insert sequence (def off)",default=False,action="store_true")
(options, args) = parser.parse_args()

def mfreq(datadict):
  mvalue = max(datadict.values())
  for key,value in datadict.items():
    if mvalue == value: return key
  return "NA"

if (options.fraction <= 0.5):
  sys.stderr.write("Error: Fraction has to be above 0.5. Exiting...")
  sys.exit()

cbarcode = None
cassignments = defaultdict(int)
cquality = defaultdict(lambda:defaultdict(int))

for line in sys.stdin:
  fields = line.rstrip().split("\t")
  if len(fields) == 3:
    tag = fields[0]
    if tag != cbarcode:
      total = float(sum(cassignments.values()))
      if (total >= options.minimum):
        reported = False
        maxcount = 0
        for insert,count in cassignments.items():
          if count > maxcount: maxcount = count
          else: continue
          if (count >= options.minimum) and (count/total >= options.fraction): 
            reported = True
            if ((insert != "other") or (options.other)):
              sys.stdout.write("%s\t%s\t%s\t%d/%d\n"%(cbarcode,insert,mfreq(cquality[insert]),count,total))
        if options.ambiguous and (not reported):
          sys.stdout.write("%s\t%s\t%s\t%d/%d\n"%(cbarcode,"ambiguous",mfreq(cquality[insert]),maxcount,total))
      if (cbarcode != None) and (tag < cbarcode):
        sys.stderr.write("Error: File is not sorted by barcode! Exiting...")
        sys.exit()
      cbarcode=tag
      cassignments = defaultdict(int)
      cquality=defaultdict(lambda:defaultdict(int))
    cassignments[fields[1]]+=1
    cquality[fields[1]][fields[2]]+=1

total = float(sum(cassignments.values()))
if (total >= options.minimum):
  reported = False
  maxcount = 0
  for insert,count in cassignments.items():
    if count > maxcount: maxcount = count
    else: continue
    if (count >= options.minimum) and (count/total >= options.fraction): 
      reported = True
      if ((insert != "other") or (options.other)):
        sys.stdout.write("%s\t%s\t%s\t%d/%d\n"%(cbarcode,insert,mfreq(cquality[insert]),count,total))
    if options.ambiguous and (not reported):
      sys.stdout.write("%s\t%s\t%s\t%d/%d\n"%(cbarcode,"ambiguous",mfreq(cquality[insert]),maxcount,total))

