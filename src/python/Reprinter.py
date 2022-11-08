#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 17:34:06 2019

@author: oshah
"""

import re, sys, time
from tabulate import tabulate
from decimal import Decimal

class Reprinter:
    def __init__(self):
        self.text = ''

    def moveup(self, lines):
        for _ in range(lines):
            sys.stdout.write("\x1b[A")

    def reprint(self, text):
        # Clear previous text by overwritig non-spaces with spaces
        self.moveup(self.text.count("\n"))
        sys.stdout.write(re.sub(r"[^\s]", " ", self.text))

        # Print new text
        lines = min(self.text.count("\n"), text.count("\n"))
        self.moveup(lines)
        sys.stdout.write(text)
        self.text = text

reprinter = Reprinter()



digits = "{:.5E}"
dat = [[digits.format(Decimal(str(0))), 
        digits.format(Decimal(str(0))),
        digits.format(Decimal(str(0))), 
        digits.format(Decimal(str(0))), 
        digits.format(Decimal(str(0)))]]

tabl = tabulate(dat, headers=["R [R_e]", "P [MPa]",\
                    'm [M_e]', 'T [K]','rho [kg m-3]'])

print (tabl)
for i in range(10):
    sys.stdout.write('\033[3A')
    dat = [[digits.format(Decimal(str(1*i))), 
        digits.format(Decimal(str(2*i))),
        digits.format(Decimal(str(3*i))), 
        digits.format(Decimal(str(4*i))), 
        digits.format(Decimal(str(5*i)))]]
    tabl = tabulate(dat, headers=["R [R_e]", "P [MPa]",\
                    'm [M_e]', 'T [K]','rho [kg m-3]'])
    print (tabl)
    time.sleep(.2)

    
print ('you suck')