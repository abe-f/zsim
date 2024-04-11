import h5py
import numpy as np
from readable_number import ReadableNumber
rn = ReadableNumber

f = h5py.File("../zsim-ev.h5", 'r')

dset = f["stats"]["root"]

endPhase = dset[-1]['phase']

totalInstrs = np.sum(dset[-1]['beefy']['instrs'])
# to get core N, do dset[-1]['beefy'][N]['instrs']
print("total instructions = " + str(rn(totalInstrs)))

l1i_hits = np.sum(dset[-1]['l1i']['hGETS'] + dset[-1]['l1i']['hGETX'] + \
                  dset[-1]['l1i']['fhGETS'] + dset[-1]['l1i']['fhGETX'] )
print("total l1i hits = " + str(rn(l1i_hits)))

l1d_hits = np.sum(dset[-1]['l1d']['hGETS'] + dset[-1]['l1d']['hGETX'] + \
                  dset[-1]['l1d']['fhGETS'] + dset[-1]['l1d']['fhGETX'] )
print("total l1d hits = " + str(rn(l1d_hits)))

print()

print("total l1 hits = " + str(rn(l1i_hits + l1d_hits)))

l2_hits = np.sum(dset[-1]['l2']['hGETS'] + dset[-1]['l2']['hGETX'])
print("total l2 hits = " + str(rn(l2_hits)))

l3_hits = np.sum(dset[-1]['l3']['hGETS'] + dset[-1]['l3']['hGETX'])
print("total l3 hits = " + str(rn(l3_hits)))

print("total cache hits = " + str(rn(l1i_hits+l1d_hits+l2_hits+l3_hits)))

# excluding mGETXSM S->M upgrade misses because those aren't really misses
l3_misses = np.sum(dset[-1]['l3']['mGETS'] + dset[-1]['l3']['mGETXIM'])
print()

print("total l3 misses = " + str(rn(l3_misses)))
print("l3 hit rate = " + str(rn(l3_hits/(l3_misses+l3_hits))))