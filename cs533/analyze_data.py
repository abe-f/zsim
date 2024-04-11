import h5py

f = h5py.File("../zsim-ev.h5", 'r')

print(list(f.keys()))

dset = f['stats']

print(dset[0])