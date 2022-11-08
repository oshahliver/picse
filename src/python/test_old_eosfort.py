import eosfort
import time

eosfort.initiate(materials=[0], ntables=1)

t0=time.time()

eosfort.interpolate(x=340.0, y=1.23e6, params=[2,3,4], ll=0)

t=time.time()

print ('elapsed time in old eosfort =', (t-t0)*1.0e3,'ms')
