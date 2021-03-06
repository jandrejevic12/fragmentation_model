import numpy as np
import os
import hashlib, hmac, pickle

def load_data(filename):
    digest, pickled = pickle.load(open(filename+'.p','rb'))
    key = bytes(filename, 'utf-8')
    # confirm data integrity
    new_digest = hmac.new(key, pickled, hashlib.sha256).hexdigest()
    if digest == new_digest:
        print("Data signature matches; data loaded.")
        return pickle.loads(pickled)
    else:
        print("Data signature does not match; data not loaded.")

def setup_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

# constant parameters
N = 3000 # image size
dim = 2830. # actual sheet size
footprint = 4
min_size = 100

ts = [[1,2,3,24],[1,2,3,24],[1,2,3,3],[1,2,3,24],[1,2,3,3],[1,2,3,3],[1,2,3,3]]
Exps = [45, 43, 37, 41, 29, 28, 39]
Deltas = [0.63,0.45,0.36,0.27,0.18,0.09,0.045]

data_prefix = 'raw_data/'
matdir = 'full_win24_20_10'
imdir = 'vector_images'

# histogram plotting
a = 1.55
bins = a**np.arange(-100,30,1)
ds = np.diff(bins)
sl = (bins[:-1]+bins[1:])*0.5

# for length:
a2 = 1.3
bins2 = a2**np.arange(-100,50,1)
ds2 = np.diff(bins2)
sl2 = (bins2[:-1]+bins2[1:])*0.5

# font sizes
major = 36
minor = 28

# model parameters
R = 0.135
A = 2.12271554
tau = 24.04077765964373
alpha = 9.031066751899312
c1 = 52.
c2 = 0.1
