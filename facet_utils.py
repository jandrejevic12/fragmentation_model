from math_imports import *
from image_imports import *
from plot_imports import *
from setup_imports import *

seed = 0
cmap = get_cmap(matplotlib.cm.get_cmap('twilight'), 0.1, 0.9)

def inflate(labels, border, maxiter=24):
    for iteration in range(maxiter):
        labels_inflated = ndimage.maximum_filter(labels, size=3).astype(np.int)
        for lab in np.unique(labels)[1:].astype(np.int):
            labels_inflated[labels==lab] = lab
        labels_inflated[border.astype(np.bool)] = 0
        labels = labels_inflated
    return labels

def get_region_props(labels):
    labels = labels.astype(np.float)
    lens = []
    for col in range(N):
        labs, cts = np.unique(labels[:,col], return_counts=True)
        lens += list(cts[1:])
    lens = np.array(lens)
    regions = regionprops(labels.astype(np.int))
    areas = np.zeros(len(regions))
    per = np.zeros(len(regions))
    inds = np.zeros((len(regions),2))
    for i,region in enumerate(regions):
        areas[i] = region.area
        per[i] = region.perimeter
        inds[i,:] = region.centroid
        #labels[labels==region.label] = np.log(region.area)
    return labels, areas, per, lens, inds

def shuffle_labels(labels, seed):
    np.random.seed(seed)
    cut = 2
    shuffled = np.unique(labels)[cut:].astype(np.int)
    np.random.shuffle(shuffled)
    new_labels = np.zeros(labels.shape)
    for i in range(cut):
        new_labels[labels == i] = i
    for i,lab in enumerate(np.unique(labels)[cut:].astype(np.int)):
        new_labels[labels == lab] = shuffled[i]
    return new_labels

def plot_map(raw_data, labels, border, inds, exp, t, seed, suffix):
    fig = Figure(figsize=(10,10))
    canvas = FigureCanvas(fig)
    ax = fig.add_subplot(111)
    ax.imshow(raw_data)
    # random color
    new_labels = shuffle_labels(labels, seed)
    ax.imshow(new_labels, cmap=cmap, alpha=0.8, vmin=1)
    ax.set_yticks([],''); ax.set_xticks([],'')
    ax.axis('off')
    fig.savefig("labeled_data/Exp"+str(exp)+"/im"+str(t)+"_"+suffix+".png", dpi=150, bbox_inches='tight', pad_inches=0)

def gen_facets(Exps, ts):
    Areas = {}
    Per = {}
    Lens = {}
    for j,exp in enumerate(Exps):
        Areas[exp] = []
        Per[exp] = []
        Lens[exp] = []
        for t in ts[j]:
            # read in the raw data and segmented facets
            raw_data = imread("labeled_data/Exp"+str(exp)+"/im"+str(t)+".png")
            border = imread("labeled_data/Exp"+str(exp)+"/im"+str(t)+"_border.png", as_grey=True)
            border = np.logical_not(border.astype(np.int)).astype(np.int)
            if t == 1 and exp == 28:
                f = "labeled_data/Exp"+str(exp)+"/D"+str(t)+"_shifted.png"
                border[:-28,:] = border[28:,:]
                raw_data[:-28,:,:] = raw_data[28:,:,:]
            elif t ==1 and exp == 39:
                f = "labeled_data/Exp"+str(exp)+"/D"+str(t)+"_shifted.png"
                border[:-7,:] = border[7:,:]
                raw_data[:-7,:,:] = raw_data[7:,:,:]
            else:
                f = "labeled_data/Exp"+str(exp)+"/D"+str(t)+".png"
            im = imread(f, as_grey=True)

            # label individual facets
            labels = label(im.astype(np.int), connectivity=2)
            
            # inflate labeled regions to eliminate crease width
            labels = inflate(labels, border)
            print(exp,t,'filter complete!')

            # get region properties
            labels, areas, per, lens, inds = get_region_props(labels)
            Areas[exp] += [areas]
            Per[exp] += [per]
            Lens[exp] += [lens]

            # plot
            plot_map(raw_data, labels, border, inds, exp, t, seed, 'manual')

    # save properties
    pickle_file = open('facet.p', 'wb')
    pickle.dump([Areas, Per, Lens], pickle_file)
    pickle_file.close()

def gen_facets_watershed(matdir, Exps, ts):
    Areas = {}
    Per = {}
    Lens = {}
    for j,exp in enumerate(Exps):
        Areas[exp] = []
        Per[exp] = []
        Lens[exp] = []
        for t in ts[j]:
            # read in the raw data and clean skeletonization
            raw_data = imread("labeled_data/Exp"+str(exp)+"/im"+str(t)+".png")
            border = imread("labeled_data/Exp"+str(exp)+"/im"+str(t)+"_border.png", as_grey=True)
            border = np.logical_not(border.astype(np.int)).astype(np.int)
            mat_contents = sio.loadmat("/scratch/shared/"+matdir+"/exp"+str(exp)+"/C"+str(t)+".mat")
            pos, neg = mat_contents["C"+str(t)+"_pos"], mat_contents["C"+str(t)+"_neg"]
            im = (pos+np.abs(neg)).astype(np.int)

            # dilate the skeletonization before labeling
            for dil in range(4):
                im = binary_dilation(im)
            im = ((im>0)|border.astype(np.bool)).astype(np.int)

            # label individual facets
            labels, markers = run_watershed(im.astype(np.int), border, footprint)
            
            # inflate labeled regions to eliminate crease width
            labels = inflate(labels, border)
            print(exp,t,'filter complete!')

            # get region properties
            labels, areas, per, lens, inds = get_region_props(labels)
            Areas[exp] += [areas]
            Per[exp] += [per]
            Lens[exp] += [lens]

            # plot
            plot_map(raw_data, labels, border, inds, exp, t, seed, 'watershed')
    
    # save properties
    pickle_file = open('facet_watershed.p', 'wb')
    pickle.dump([Areas, Per, Lens], pickle_file)
    pickle_file.close()

def run_watershed(im, border, footprint):
    N = im.shape[1]
    image = np.logical_not(im)
    distance = ndimage.distance_transform_edt(image)
    local_maxi = h_maxima(distance, footprint)
    local_maxi[0,:] = 0; local_maxi[-1,:] = 0; local_maxi[:,0] = 0; local_maxi[:,-1] = 0;
    markers = ndimage.label(local_maxi)[0]
    labels = watershed(-distance, markers, mask=image)
    labels = labels.astype(np.float)
    for i in np.unique(labels):
        inds = np.argwhere(labels==i)
        if np.any(inds[:,0]==0) or np.any(inds[:,0]==N) or np.any(inds[:,1]==0) or np.any(inds[:,1]==N):
            labels[labels==i] = -1
    labels[labels==-1] = 0
    labels[border.astype(np.bool)] = 0
    return labels, markers

if __name__ == "__main__":
    ts = [[1,2,3,24],[1,2,3,24],[1,2,3,3],[1,2,3,24],[1,2,3,3],[1,2,3,3],[1,2,3,3]]
    # extension using automated method:
    ts = [[int(i+1) for i in range(24)],
          [int(i+1) for i in range(24)],
          [int(i+1) for i in range(24)],
          [int(i+1) for i in range(24)],
          [int(i+1) for i in range(19)],
          [int(i+1) for i in range(10)],
          [int(i+1) for i in range(4)]]
    Exps = [45, 43, 37, 41, 29, 28, 39]
    #gen_facets(Exps, ts)
    gen_facets_watershed('full_win24_20_10_ext', Exps, ts)
