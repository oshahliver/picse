import numpy as np
import fortplanet
import PlanetFort
import time
import argparse
from astropy.io import ascii
import astropy
        
parser = argparse.ArgumentParser(description='Run Population')
parser.add_argument('-N', type=int, default=10)

args = parser.parse_args()
args_dict = vars(args)

#fortplanet.wrapper.load_eos_tables(table_dir='/home/os18o068/Documents/PHD/Projects/Planets/eos_tables/')

t0=time.time()

fortplanet.wrapper.do_stuff(n_planets=args_dict['N'])

t = time.time()

print ("Elapsed time in fortplanet:", round((t-t0)*1.0e3,2), "ms")

