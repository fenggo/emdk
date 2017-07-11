#!/usr/bin/env python
# -*- coding: utf-8 -*-
from os import system
from os.path import isfile
from commands import getoutput

fs = open('structureDatabase.f90','r')
fsr= open('structureDatabase.r.f90','w')

lread = False
for line in fs.readlines():                # prepare input file
    #if len(line.split())>=1:
       if line.find('types = (/ &')>=0:
          lread = True
          nt = 1
       elif line.find('/)')>=0:
          lread = False  
       if not lread and line.find('  /)')!=11:
          print>>fsr,'%s' %line[:-1] # delete /n                           
       elif lread and line.find('types')<0:
          ll = line.split()[0]
          if len(ll)>1:
             ll = ll[:-1]
          print>>fsr, '    types(',nt,') = ',ll
          nt += 1

fs.close()
fsr.close()

