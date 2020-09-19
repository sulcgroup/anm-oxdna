import Bio
import Bio.PDB
import scipy.linalg
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as la
import math
import sys
import os

import copy


#Set cuda_support to False if CUDA if not needed
#To enable cuda support must install pycuda and skcuda
#AND have LD_LIBRARY_PATH variable set correctly
#might as well set PATH while you're at it

# Ubuntu Example
# ex. PATH=/usr/local/cuda/bin
# ex. LD_LIBRARY_PATH=/usr/local/cuda/lib64
cuda_support = False #Enable CUDA operations by loading CUDA drivers and initializing GPU

if cuda_support:
    from pycuda import autoinit
    from pycuda import gpuarray
    import pycuda.driver as cuda
    import skcuda.linalg as cuda_la
    import skcuda.misc as cuda_misc
    cuda_misc.init()
    cuda_la.init()


# Helper Functions
def dist(coord, i, j):
    dx = coord[i, 0] - coord[j, 0]
    dy = coord[i, 1] - coord[j, 1]
    dz = coord[i, 2] - coord[j, 2]
    return math.sqrt(dx ** 2 + dy ** 2 + dz ** 2)

def flatten(l1):
    return [item for sublist in l1 for item in sublist]

def divide_list(list, piecenum):
    divis = len(list) % piecenum
    if divis == 0:
        diff = len(list) // piecenum
    else:
        print("List is not divisble by the number of pieces you provided")
        sys.exit()

    for i in range(piecenum):
        if i == piecenum - 1:
            sublist = list[i * diff:]
            yield (sublist)
        else:
            sublist = list[i * diff:(i + 1) * diff]
            yield (sublist)

def normalize_vector(v):
    if v[0] == 0. and v[1] == 0. and v[2] == 0.:
        return v
    else:
        return np.divide(v, np.sqrt(np.sum(v**2)))

def map_bfactors(adjseq, exp_bfacts, fillvalue=1.):
    bfmap = []
    count = 0
    for xid, x in enumerate(list(adjseq)):
        if x == '-':
            count += 1
            bfmap.append(fillvalue)
        else:
            bfmap.append(exp_bfacts[xid - count])
    return bfmap


def read_seqfile(seqfile):
    o = open(seqfile, 'r')
    data = o.read()
    o.close()
    data.split('\n')
    datastr = ''.join(data)
    newdata = [x for x in list(datastr) if x.isalpha() or x =='-']
    fullseq = ''.join(newdata)
    return fullseq



conv = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F',
        'ASN': 'N',
        'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y',
        'MET': 'M'}

# Gets all chains in PDB File
# Return Types

def get_pdb_info(pdb_file, returntype='cb'):
    print("INFO: Don't worry about Discontinuous chain warnings!")
    if "/" in pdb_file:
        pdbid = pdb_file.rsplit('/', 1)[1].split('.')[0]
    else:
        pdbid = pdb_file.split('.')[0]
    structure = Bio.PDB.PDBParser().get_structure(pdbid, pdb_file)

    #Get all chains in file
    model = Bio.PDB.Selection.unfold_entities(structure, 'C')

    #orientation vector container
    tmp_N_vectors = []

    chainids, chain_coords, chain_seqs, chain_bfactors, chainbounds = [], [], [], [], {}

    # iterate through chains in pdb
    for chain in model:
        # Add chain
        chainids.append(chain.get_id())
        # empty coordinates each iteration to separate by chain
        coordtmp = []
        bfacttmp = []
        chainseqtmp = []
        chainboundstmp = []
        # iterate through residues
        for residue in chain.get_residues():
            tags = residue.get_full_id()
            # print(tags)
            if tags[3][0] == " ":
                # Get Residues one letter code
                onelettercode = conv[residue.get_resname()]
                # get residue number and identity per chain
                chainseqtmp.append((tags[2], onelettercode))
                chainboundstmp.append(int(tags[3][1]))
                atoms = residue.get_atoms()
                # Center of Mass Used only if Orientation=True
                com = np.full(3, 0.0)
                count = 0
                for atom in atoms:
                    if atom.get_id() == 'CA':
                        coordinates = atom.get_coord()
                        bfactor = atom.get_bfactor()
                        bfacttmp.append(bfactor)
                        # Add coordinates to tmp container
                        coordtmp.append(coordinates)
                        if atom.get_id() != 'N' or atom.get_id() != 'C' or atom.get_id() != 'O':
                            count += 1
                            com += atom.get_coord()
                if np.sum(com**2) != 0.:
                    com /= count
                    # Vector from CA coordinates to Center of Mass
                    nvec = normalize_vector(com - coordtmp[-1])
                    tmp_N_vectors.append(nvec)
                else:
                    tmp_N_vectors.append(np.asarray([0., 0., 0.]))
        # Before next chain add all of this chain's coordinates to chain_coords
        chainbounds[str(tags[2])] = (min(chainboundstmp), max(chainboundstmp))
        chain_seqs.append(chainseqtmp)
        chain_coords.append(coordtmp)
        chain_bfactors.append(bfacttmp)

    chainmap = []  # Stores missing res. and chain takeover points
    c = 0.
    for chain in chain_seqs:
        tag = chain[0][0]
        indx = c
        chainmap.append((tag, indx))
        c += len(chain)

    if(len(chain_coords) > 1):
        cflat = flatten(chain_coords)
        N = len(cflat)
        print('INFO: N =', N)

    # Dealing with Null Normal Vectors
    normal_vectors = []
    fcoord = np.asarray(flatten(chain_coords))
    for iid, i in enumerate(tmp_N_vectors):
        if np.sum(i**2) == 0. and iid != len(tmp_N_vectors)-1 and iid != 0:
            #vector from particle i to i + 1
            rij = fcoord[iid+1] - fcoord[iid]
            #vector from particle i to i-1
            rih = fcoord[iid-1] - fcoord[iid]
            nvec = normalize_vector(np.cross(rij, rih))
            normal_vectors.append(nvec)
        elif np.sum(i**2) == 0:
            nvec = np.asarray([1., 0., 0.])
            normal_vectors.append(nvec)
        else:
            nvec = normalize_vector(i)
            normal_vectors.append(nvec)
    a1s = np.asarray(normal_vectors)

    # Calc a3 vectors
    basis_vecs = np.asarray([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
    # Threshold to declare whether vector is perpendicular or parallel
    epsilon = 1e-3
    tmp_orth_vectors = []
    for a1 in a1s:
        if np.linalg.norm(a1) != 1.:
            a1 = normalize_vector(a1)
        # check if parallel to 1st basis vector
        if not (np.dot(a1, basis_vecs[0]) > 1 - epsilon):
            a3 = np.cross(a1, basis_vecs[0])
        elif not (np.dot(a1, basis_vecs[1]) > 1 - epsilon):
            a3 = np.cross(a1, basis_vecs[1])
        # Arbitrarily Defined Just needs to be perpendicular
        tmp_orth_vectors.append(a3)
    a3s = np.asarray(tmp_orth_vectors)

    fcoords = np.asarray(flatten(chain_coords))
    fbfact = np.asarray(flatten(chain_bfactors))
    fseq = []
    for i in chain_seqs:
        fseq += i

    returndict = {"b":fbfact, "B":chain_bfactors, "c":fcoords, "C":chain_coords, "s":chain_seqs, "S": fseq,
                  "m":chainmap, "1":a1s, "3":a3s, "p":chainbounds}

    returnstrlist = list(returntype)
    returnparser = [returndict[i] for i in returnstrlist]
    if len(returnparser) > 1:
        return tuple(returnparser)
    else:
        return returnparser[0]


def free_compare(out, *data, legends=[], title='', bmap=''):
    cs = ['r', 'c', 'b', 'm', 'k', 'y']
    ls = ['-', '--', '-.', ':']
    fig, ax = plt.subplots(1)
    ax.set_title(title)
    plt.xlabel('Residue', fontsize=16)
    plt.ylabel(r'Mean Squared Deviation [$\AA^2$]', fontsize=16)
    for yid, y in enumerate(data):
        if legends:
            ax.plot(np.arange(0, len(y), 1), y, c=cs[yid % 6], label=legends[yid], linestyle=ls[yid % 4])
        else:
            ax.plot(np.arange(0, len(y), 1), y, c=cs[yid % 6], linestyle=ls[yid % 4])
    if bmap:
        chainids, indxs = zip(*bmap)
        for i in range(len(chainids)):
            plt.axvline(indxs[i])
            plt.text(indxs[i]+1,0, chainids[i].upper(), rotation=90, fontsize=14)
    ax.legend(loc=2, prop={'size': 14})
    plt.savefig(out, dpi=600)
    plt.close()

def free_compare_jupyter(*data, legends=[], title='', bmap=''):
    cs = ['r', 'c', 'b', 'm', 'k', 'y']
    ls = ['-', '--', '-.', ':']
    fig, ax = plt.subplots(1)
    ax.set_title(title)
    plt.xlabel('Residue', fontsize=16)
    plt.ylabel(r'Mean Squared Deviation [$\AA^2$]', fontsize=16)
    for yid, y in enumerate(data):
        if legends:
            ax.plot(np.arange(0, len(y), 1), y, c=cs[yid % 6], label=legends[yid], linestyle=ls[yid % 4])
        else:
            ax.plot(np.arange(0, len(y), 1), y, c=cs[yid % 6], linestyle=ls[yid % 4])
    if bmap:
        chainids, indxs = zip(*bmap)
        for i in range(len(chainids)):
            plt.axvline(indxs[i])
            plt.text(indxs[i]+1,0, chainids[i].upper(), rotation=90, fontsize=14)
    ax.legend(loc=2, prop={'size': 14})



def diff_sqrd(l1, l2):
    diff = 0
    for i in range(len(l1)):
        diff += abs(l1[i] ** 2 - l2[i] ** 2)
    return diff

def load_sim_rmsds_from_file(file):
    o = open(file, 'r')
    str_arr = o.read()
    o.close()
    # Removes beginning -> '{"RMSD (nm)"":
    arr1 = (str_arr.split(':')[1]).lstrip()
    # Removes End
    arr2 = arr1.split('}')[0]
    # Convert to List of Strings
    arr3 = arr2.strip('][]').split(', ')
    # Convert to List of Floats
    rmsds = [float(x) for x in arr3]
    sim_msds = [x ** 2 * 100 for x in rmsds]
    sim_bfactors = [(8 * math.pi ** 2) / 3 * x ** 2 * 100 for x in rmsds]
    return sim_msds, sim_bfactors

def get_spherical_angles(v1):
    if np.sqrt(np.sum(v1**2)) != 1:
        v0 = normalize_vector(v1)
    else:
        v0 = v1
    #In radians
    phi = math.atan2(v0[1], v0[0])
    theta = math.acos(v0[2])
    return theta, phi

def spherical_dtheta(theta, phi):
    return np.asarray([math.cos(theta)*math.cos(phi), math.cos(theta)*math.sin(phi), -1.*math.sin(theta)])

def spherical_dphi(theta, phi):
    return np.asarray([-1.*math.sin(theta)*math.sin(phi), math.cos(phi)*math.sin(theta), 0.])


class ANM(object):

    def __init__(self, coord, exp_bfactors, T=300, cutoff=15):
        #Conversion factors for B Factors to Fluctiations
        self.bconv = (8. * math.pi ** 2) / 3.
        self.ibconv = 3./(8. * math.pi ** 2)
        #These are important for most loops
        self.sim_length = 8.518 # A/sim length
        self.coord = np.asarray(coord)
        self.cc = len(coord)
        #Data to match to
        self.exp_bfactors = exp_bfactors
        self.exp_msds = [self.ibconv*x for x in self.exp_bfactors]

        self.msds = []
        self.ana_bfactors = []
        self.ana_msd = []
        self.ana_gamma = 0.
        # Angstroms
        self.cutoff = cutoff
        #Angstroms in 1 sim unit length
        self.sim_force_const = .05709  # (Sim Units to pN/A)
        # IN picoNetwtons/ Angstroms
        self.kb = 0.00138064852
        #Kelvin
        self.T = T

        self.distance_matrix = []
        self.model_id = 'ANM'

    def calc_dist_matrix(self):
        d_matrix = np.full((self.cc, 4*self.cc), 0.0, dtype=np.float32)
        for i in range(self.cc):
            for j in range(self.cc):
                if i >= j:
                    continue
                else:
                    #Calculating distance
                    dx = self.coord[j,0]- self.coord[i,0]
                    dy = self.coord[j,1]- self.coord[i,1]
                    dz = self.coord[j,2]- self.coord[i,2]
                    # dx = self.coord[i, 0] - self.coord[j, 0]
                    # dy = self.coord[i, 1] - self.coord[j, 1]
                    # dz = self.coord[i, 2] - self.coord[j, 2]
                    dist = np.sqrt(dx**2 + dy**2 + dz**2)
                    # too far, skips
                    if dist > float(self.cutoff):
                        continue
                    else:
                        d_matrix[i, 4*j+1] = dx
                        d_matrix[i, 4*j+2] = dy
                        d_matrix[i, 4*j+3] = dz
                        d_matrix[i, 4*j] = dist
        self.distance_matrix = d_matrix

    def calc_hess_fast_sc(self, spring_constant_matrix):
        threeN = 3 * self.cc
        hess = np.zeros((threeN, threeN), dtype=np.float32)
        for i in range(self.cc):
            for j in range(self.cc):
                if i >= j:
                    continue
                di = self.distance_matrix[i, 4 * j:4 * j + 4]
                # Filter so that Hessian is only created for those bonds in bonds array
                if di[0] != 0:
                    di2 = np.square(di)
                    g = spring_constant_matrix[i, j]

                    diag = g * di2[1:4] / di2[0]

                    xy = g * (di[1] * di[2]) / di2[0]
                    xz = g * (di[1] * di[3]) / di2[0]
                    yz = g * (di[2] * di[3]) / di2[0]

                    full = np.asarray([[diag[0], xy, xz], [xy, diag[1], yz], [xz, yz, diag[2]]], order='F')

                    # Hii and Hjj
                    hess[3 * i:3 * i + 3, 3 * i:3 * i + 3] += full
                    hess[3 * j: 3 * j + 3, 3 * j:3 * j + 3] += full

                    # Hij and Hji
                    hess[3 * i: 3 * i + 3, 3 * j: 3 * j + 3] -= full
                    hess[3 * j: 3 * j + 3, 3 * i: 3 * i + 3] -= full
        return hess

    def calc_hess_fast_unitary(self, gamma=1.):
        threeN = 3 * self.cc
        hess = np.zeros((threeN, threeN), dtype=np.float32)
        for i in range(self.cc):
            for j in range(self.cc):
                if i >= j:
                    continue
                di = self.distance_matrix[i, 4 * j:4 * j + 4]
                # Filter so that Hessian is only created for those bonds in bonds array
                if di[0] != 0:
                    di2 = np.square(di)
                    g = gamma

                    diag = g * di2[1:4] / di2[0]

                    xy = g * (di[1] * di[2]) / di2[0]
                    xz = g * (di[1] * di[3]) / di2[0]
                    yz = g * (di[2] * di[3]) / di2[0]

                    full = np.asarray([[diag[0], xy, xz], [xy, diag[1], yz], [xz, yz, diag[2]]], order='F')

                    # Hii and Hjj
                    hess[3 * i:3 * i + 3, 3 * i:3 * i + 3] += full
                    hess[3 * j: 3 * j + 3, 3 * j:3 * j + 3] += full

                    # Hij and Hji
                    hess[3 * i: 3 * i + 3, 3 * j: 3 * j + 3] -= full
                    hess[3 * j: 3 * j + 3, 3 * i: 3 * i + 3] -= full
        return hess

    def calc_inv_Hess(self, hess, cuda=False):
        if cuda:
            try:
                thess = np.asarray(hess, dtype=np.float32, order='C')
                cu_hess = gpuarray.to_gpu(thess)
                u_gpu, s_gpu, vh_gpu = cuda_la.svd(cu_hess, 'S', 'S')
                U, w, Vt = u_gpu.get(), s_gpu.get(), vh_gpu.get()
            except NameError:
                print('CUDA failed to Initialize properly. Using CPU version')
                U, w, Vt = scipy.linalg.svd(hess, full_matrices=False)
        else:
            U, w, Vt = scipy.linalg.svd(hess, full_matrices=False)
        #print(w)
        S = scipy.linalg.diagsvd(w, len(w), len(w))
        tol = 1e-4
        singular = w < tol
        invw = 1 / w
        invw[singular] = 0.
        hessinv = np.dot(np.dot(U, np.diag(invw)), Vt)
        return hessinv

    def save_inverse_Hessian(self, invhess, outfile):
        np.save(outfile, invhess)

    def load_inverse_Hessian(self, infile):
        invhess = np.load(infile)
        return invhess

    #finds closest particles to the coordindx
    # def find_nearest_particles(self, coordindx, cutoff):
    #
    #
    # def anm_remove_peak(self):
    #     hpeak = np.argmax(np.asarray(self.msds))

    def ANM_search(self, br, er, step):
        r = []
        sc_range = np.arange(br, er, step)
        for i in sc_range:
            g = float(i)
            r.append(diff_sqrd([x * 1. / g * self.bconv for x in self.msds], self.exp_bfactors))

        results = np.array(r)
        bg = np.argmin(results)
        return bg * step + br

    def ANM_fit_to_exp(self, start=0.001, end=5.000, step=0.001):
        g1 = self.ANM_search(start, end, step)
        self.ana_gamma = g1
        self.ana_msd = [x * 1 / self.ana_gamma for x in self.msds]
        self.ana_bfactors = [self.bconv * x * 1 / self.ana_gamma for x in self.msds]

    def ANM_fit_to_exp_linear(self):
        if self.msds: # The mean square deviations must be precomputed
            try:
                from sklearn.linear_model import LinearRegression
            except ImportError:
                print('Check that sklearn module is installed')
                sys.exit()
            flex_data = np.asarray([x * self.bconv for x in self.msds])
            exp_data = np.asarray(self.exp_bfactors)
            X = flex_data.reshape(-1, 1) #.transpose()
            Y = exp_data
            fitting = LinearRegression(fit_intercept=False)
            fitting.fit(X, Y)
            slope = fitting.coef_
            #print(slope)
            #print(X, Y)
            self.ana_gamma = float(1/slope)
            self.ana_msd = [x * 1 / self.ana_gamma for x in self.msds]
            self.ana_bfactors = [self.bconv * x * 1 / self.ana_gamma for x in self.msds]
        else:
            print("Prior to Calling ANM_fit_to_exp_linear you must calc the msds")

    def calc_msds(self, invhess):
        self.msds = []
        for i in range(self.cc):
            self.msds.append(self.kb * self.T * (invhess[3 * i, 3 * i] + invhess[3 * i + 1, 3 * i + 1] +
                                                 invhess[3 * i + 2, 3 * i + 2]))

    def anm_calc_bfactors(self):
        self.ana_bfactors = []
        self.ana_bfactors = [self.bconv*x*1/self.ana_gamma for x in self.ana_msd]

    def calc_ANM_sc(self, spring_constant_matrix, cuda=False):
        self.calc_dist_matrix()
        hess = self.calc_hess_fast_sc(spring_constant_matrix)
        iH = self.calc_inv_Hess(hess, cuda=cuda)
        self.calc_msds(iH)
        self.ana_bfactors = [self.bconv * x for x in self.msds]

    def calc_ANM_unitary(self, cuda=False):
        self.calc_dist_matrix()
        hess = self.calc_hess_fast_unitary()
        iH = self.calc_inv_Hess(hess, cuda=cuda)

        self.calc_msds(iH)
        self.ANM_fit_to_exp_linear()

    def anm_compare_bfactors(self, outfile, bmap=''):
        if self.ana_bfactors:
            if bmap:
                free_compare(outfile, self.exp_bfactors, self.ana_bfactors, bmap=bmap,
                             legends=['Experimental  (PDB)',
                                      'Analytical (ANM)' + str(round(self.ana_gamma * 100, 3)) + "(pN/A)"])
            else:
                free_compare(outfile, self.exp_bfactors, self.ana_bfactors,
                    legends=['Experimental  (PDB)', 'Analytical (ANM)' + str(round(self.ana_gamma*100, 3))+ "(pN/A)"])
        else:
            print('Analytical B Factors have not been Calculated')

    def anm_compare_bfactors_jupyter(self, bmap=''):
        if self.ana_bfactors:
            if bmap:
                free_compare_jupyter(self.exp_bfactors, self.ana_bfactors, bmap=bmap,
                             legends=['Experimental  (PDB)',
                                      'Analytical (ANM)' + str(round(self.ana_gamma * 100, 3)) + "(pN/A)"])
            else:
                free_compare_jupyter(self.exp_bfactors, self.ana_bfactors,
                             legends=['Experimental  (PDB)',
                                      'Analytical (ANM)' + str(round(self.ana_gamma * 100, 3)) + "(pN/A)"])
        else:
            print('Analytical B Factors have not been Calculated')

    def simplifty_2d_matrix(self, matrix, tol=0.001):
        # print(np.shape(matrix))
        matrix[matrix < (tol/self.sim_force_const)] = 0.
        # newmatrix = copy.deepcopy(matrix)
        # newmatrix[low_val_indices] = 0.
        return matrix


class peptide():
    def __init__(self, pdbfile, cutoff=5, indx_cutoff=5, potential=5, backbone_weight=0):
        self.model_id = 'pep'
        ccoord, cbfacts = get_pdb_info(pdbfile, returntype='cb')
        self.coord = ccoord
        self.pi = len(self.coord)
        self.bonds = []
        self.cutoff = cutoff
        self.indx_cutoff = indx_cutoff
        self.potential = potential
        self.backbone_weight=backbone_weight
    def calc_bonds_bounded(self, indxLowerbound, indxUpperbound):
        self.bonds = []
        for i in range(self.pi):
            if i < self.pi-1:
                self.bonds.append((i, i+1))
            for j in range(self.pi):
                if i >= j:
                    continue
                if i < indxLowerbound or j < indxLowerbound:
                    continue
                if i > indxUpperbound or j > indxUpperbound:
                    continue
                dx = self.coord[j, 0] - self.coord[i, 0]
                dy = self.coord[j, 1] - self.coord[i, 1]
                dz = self.coord[j, 2] - self.coord[i, 2]
                dist = math.sqrt(abs(dx) ** 2 + abs(dy) ** 2 + abs(dz) ** 2)
                if dist < self.cutoff:
                    self.bonds.append((i, j))
    def calc_bonds(self):
        self.bonds = []
        for i in range(self.pi):
            for j in range(self.pi):
                if i >= j:
                    continue
                if j > i+self.indx_cutoff:
                    continue
                dx = self.coord[j, 0] - self.coord[i, 0]
                dy = self.coord[j, 1] - self.coord[i, 1]
                dz = self.coord[j, 2] - self.coord[i, 2]
                dist = math.sqrt(abs(dx) ** 2 + abs(dy) ** 2 + abs(dz) ** 2)
                if dist < self.cutoff:
                    self.bonds.append((i, j))


class kernel():
    def __init__(self, id):
        self.id = id
        #bfactor fitting coefficients
        self.fitting_a = 0
        # Main Storage
        self.kerns = []
        self.mu = []
        self.cutoff = 0

    def algorithim(self, dist, scale_resolution, k_factor, alg='ge'):
        # Can choose between Generalized Exponential and Generalized Lorentz Function
        def gen_exp(dist):
            return math.exp(-1 * (dist / scale_resolution) ** k_factor)

        def gen_lor(dist):
            return 1. / (1. + (dist / scale_resolution) ** k_factor)

        if alg == 'ge':
            return gen_exp(dist)
        elif alg == 'gl':
            return gen_lor(dist)

    def construct_kernel(self, coordinates, cutoff, scale_resolution, k_factor, alg='ge'):
        self.kerns = []
        self.mu = []
        self.cutoff = cutoff
        cc = len(coordinates)
        for i in range(cc):
            ker_i = 0.
            for j in range(cc):
                d = dist(coordinates, i, j)
                if cutoff > 0. and d <= cutoff:
                    ker = self.algorithim(d, scale_resolution, k_factor, alg=alg)
                elif cutoff > 0. and d > cutoff:
                    ker = 0.
                else:
                    ker = self.algorithim(d, scale_resolution, k_factor, alg=alg)
                self.kerns.append(ker)
                ker_i += ker
            self.mu.append(ker_i)

        # replace ii with sum
        for i in range(cc):
            indx = i * cc + i
            #Negative Sign Yay or Nay? Try Nay first
            self.kerns[indx] = self.mu[i]

    def shift_mu(self, shift):
        self.mu = [x * shift for x in self.mu]

    def shift_kerns(self, shift):
        self.kerns = [x * shift for x in self.kerns]


#Known Issues: Can sometimes return negative spring constants
#No python linear fit function enforces the slope values to be non-negative AND fit without an intercept
class Multiscale_ANM(ANM):
    def __init__(self, coord, exp_bfactors, T=300):
        super().__init__(coord, exp_bfactors, T=T)
        self.model_id = 'mANM'
        self.kernels = []
        self.inv_bfactors = [1./x for x in exp_bfactors]
        self.spring_constant_matrix = []

    def make_kernel(self, id, cutoff, scale_resolution, k_factor, alg='ge'):
        k = kernel(id)
        k.construct_kernel(self.coord, cutoff, scale_resolution, k_factor, alg=alg)
        self.kernels.append(k)


    def fit_kernels(self):
        raw = [k.mu for k in self.kernels]
        try:
            from sklearn.linear_model import LinearRegression
        except ImportError:
            print('Check that sklearn module is installed')
            sys.exit()
        flex_data = np.asarray(raw)
        exp_data = np.asarray(self.inv_bfactors)
        X = flex_data.transpose()
        Y = exp_data
        print(flex_data)
        fitting = LinearRegression(fit_intercept=False)
        fitting.fit(X, Y)
        slope = fitting.coef_
        print(slope)
        for i in range(len(self.kernels)):
            self.kernels[i].fitting_a = slope[i]

    def calc_msds(self, invhess):
        self.msds = []
        for i in range(self.cc):
            self.msds.append(invhess[3 * i, 3 * i] + invhess[3 * i + 1, 3 * i + 1] +
                             invhess[3 * i + 2, 3 * i + 2])

    def calc_ANM_sc(self, spring_constant_matrix, cuda=False):
        self.calc_dist_matrix()
        hess = self.calc_hess_fast_sc(spring_constant_matrix)
        iH = self.calc_inv_Hess(hess, cuda=cuda)
        self.calc_msds(iH)
        self.ana_bfactors = [x * self.ibconv for x in self.msds]


    def mANM_theor_bfactors(self, outfile, cuda=False):
        for k in self.kernels:
            print(k.kerns[0])
            k.shift_kerns(k.fitting_a)
            print(k.kerns[0])

        self.cutoff = max([k.cutoff for k in self.kernels])
        print(self.cutoff)

        spring_constant_total = np.full((self.cc, self.cc), 0.)
        for i in range(self.cc):
            for j in range(self.cc):
                if i==j:
                    k_ij = 1.
                    spring_constant_total[i, j] = k_ij
                elif i>j:
                    continue
                else:
                    k_ij = 0.
                    for k in self.kernels:
                        k_ij += k.kerns[i * self.cc + j]
                    spring_constant_total[i, j] = k_ij
                    spring_constant_total[j, i] = k_ij
        self.spring_constant_matrix = spring_constant_total

        self.calc_ANM_sc(self.spring_constant_matrix, cuda=cuda)
        self.anm_compare_bfactors(outfile)



class ANMT(ANM):
    def __init__(self, coord, exp_bfactors, T=300, cutoff=13):
        super().__init__(coord, exp_bfactors, T=T, cutoff=cutoff)
        self.model_id = 'ANMT'
        self.kernels = []

#     def get_rotation_matrix(self, axis, anglest):
#         # copied from oxDNA UTILS (Not currently used but could be helpful later)
#         # the argument anglest can be either an angle in radians
#         # (accepted types are float, int or np.float64 or np.float64)
#         # or a tuple [angle, units] where angle a number and
#         # units is a string. It tells the routine whether to use degrees,
#         # radiants (the default) or base pairs turns
#         if not isinstance(anglest, (np.float64, np.float32, float, int)):
#             if len(anglest) > 1:
#                 if anglest[1] in ["degrees", "deg", "o"]:
#                     # angle = np.deg2rad (anglest[0])
#                     angle = (np.pi / 180.) * (anglest[0])
#                 elif anglest[1] in ["bp"]:
#                     angle = int(anglest[0]) * (np.pi / 180.) * (35.9)
#                 else:
#                     angle = float(anglest[0])
#             else:
#                 angle = float(anglest[0])
#         else:
#             angle = float(anglest)  # in degrees, I think
#
#         axis = np.array(axis)
#         axis /= np.sqrt(np.dot(axis, axis))
#
#         ct = np.cos(angle)
#         st = np.sin(angle)
#         olc = 1. - ct
#         x, y, z = axis
#
#         return np.array([[olc * x * x + ct, olc * x * y - st * z, olc * x * z + st * y],
#                          [olc * x * y + st * z, olc * y * y + ct, olc * y * z - st * x],
#                          [olc * x * z - st * y, olc * y * z + st * x, olc * z * z + ct]])
#
#     def calc_bend_hess(self, kb=1, kt=1):
#         # prototype, this is not working yet
#         _7n = 7*len(self.coord)
#         self.calc_dist_matrix()
#         bend_hess = np.full((_7n, _7n), 0.0)
#         for i in range(self.cc-1):
#             j = i+1
#             a1 = self.normal_vectors[i]
#             b1 = self.normal_vectors[j]
#             a3 = self.a3s[i]
#             b3 = self.a3s[j]
#
#             r, dx, dy, dz = self.distance_matrix[i, 4 * j:4 * j + 4]
#             rij = np.asarray([dx, dy, dz])
#             vx = np.asarray([dx**2., dx*dy, dx*dz])
#             vy = np.asarray([dx*dy, dy**2., dz*dy])
#             vz = np.asarray([dx*dz, dy*dz, dz**2.])
#
#             Aix = a1[0]/r - np.dot(a1, vx)/(r**3.)
#             Aiy = a1[1]/r - np.dot(a1, vy)/(r**3.)
#             Aiz = a1[2]/r - np.dot(a1, vz)/(r**3.)
#
#             Ajx = b1[0]/r - np.dot(b1, vx)/(r**3.)
#             Ajy = b1[1]/r - np.dot(b1, vy)/(r**3.)
#             Ajz = b1[2]/r - np.dot(b1, vz)/(r**3.)
#
#             # Spherical Vector Derivatives
#             t1i, p1i = get_spherical_angles(a1)
#             t1j, p1j = get_spherical_angles(b1)
#             t3i, p3i = get_spherical_angles(a3)
#             t3j, p3j = get_spherical_angles(b3)
#
#             # Needed for a1, a3 dof derivatives
#             dv1i_theta = spherical_dtheta(t1i, p1i)
#             dv1i_phi = spherical_dphi(t1i, p1i)
#
#             dv1j_theta = spherical_dtheta(t1j, p1j)
#             dv1j_phi = spherical_dphi(t1j, p1j)
#
#             dv3i_theta = spherical_dtheta(t3i, p3i)
#             dv3i_phi = spherical_dphi(t3i, p3i)
#
#             dv3j_theta = spherical_dtheta(t3j, p3j)
#             dv3j_phi = spherical_dphi(t3j, p3j)
#
#             # Derivatives without Spring Interaction for xi, xj -> zi, zj
#             # -dir for Hij and Hji, + for Hii and Hjj
#             dir = np.full((3, 3), 0.0)
#             dir[0] = [Aix*Aix+Ajx*Ajx, Aix*Aiy+Ajx*Ajy, Aix*Aiz+Ajx*Ajz]
#             dir[1] = [Aiy*Aix+Ajy*Ajx, Aiy*Aiy+Ajy*Ajy, Aiy*Aiz+Ajy*Ajz]
#             dir[2] = [Aiz*Aix+Ajz*Ajx, Aiz*Aiy+Ajz*Ajy, Aiz*Aiz+Ajz*Ajz]
#             # Just Missing a kb term now
#
#             #Hij Angular Derivatives
#             # Angular Derivatives 3 Rows(xi, yi, zi) Two Columns(theta1j, phi1j)
#             # Postion -> Angle PA
#             #ex. dU/ dxi theta1j ...
#             PA_ij = np.full((3, 2), 0.0)
#             PA_ij[0] = [-1.*np.dot(rij, dv1j_theta)/r*Ajx, -1.*np.dot(rij, dv1j_phi)/r*Ajx]
#             PA_ij[1] = [-1.*np.dot(rij, dv1j_theta)/r*Ajy, -1.*np.dot(rij, dv1j_phi)/r*Ajy]
#             PA_ij[2] = [-1.*np.dot(rij, dv1j_theta)/r*Ajz, -1.*np.dot(rij, dv1j_phi)/r*Ajz]
#             # Just Missing a kb term now
#
#             # Angular Derivatives 2 Rows(theta1i, phi1i) Three Columns(xj, yj, zj)
#             # Angle -> Position AP
#             # ex. dU/ dtheta1i xj ...
#             AP_ij = np.full((2, 3), 0.0)
#             AP_ij[0] = [np.dot(rij, dv1i_theta) / r * Aix, np.dot(rij, dv1i_theta) / r * Aiy, np.dot(rij, dv1i_theta) / r * Aiz]
#             AP_ij[1] = [np.dot(rij, dv1i_phi) / r * Aix, np.dot(rij, dv1i_phi) / r * Aiy, np.dot(rij, dv1i_phi) / r * Aiz]
#             # Just Missing a kb term now
#
#             # Hji Angular Derivatives
#             # Angular Derivatives 3 Rows(xi, yi, zi) Two Columns(theta1j, phi1j)
#             # Postion -> Angle PA
#             # ex. dU/ dxj theta1i ...
#             PA_ji = np.full((3, 2), 0.0)
#             PA_ji[0] = [np.dot(rij, dv1i_theta) / r * Aix, np.dot(rij, dv1i_phi) / r * Aix]
#             PA_ji[1] = [np.dot(rij, dv1i_theta) / r * Aiy, np.dot(rij, dv1i_phi) / r * Aiy]
#             PA_ji[2] = [np.dot(rij, dv1i_theta) / r * Aiz, np.dot(rij, dv1i_phi) / r * Aiz]
#             # Just Missing a kb term now
#
#             # Angular Derivatives 2 Rows(theta1i, phi1i) Three Columns(xj, yj, zj)
#             # Angle -> Position AP
#             # ex. dU/ dtheta1i xi ...
#             AP_ji = np.full((2, 3), 0.0)
#             AP_ji[0] = [-1.*np.dot(rij, dv1j_theta) / r * Ajx, -1.*np.dot(rij, dv1j_theta) / r * Ajy, -1.*np.dot(rij, dv1j_theta) / r * Ajz]
#             AP_ji[1] = [-1.*np.dot(rij, dv1j_phi) / r * Ajx, -1.*np.dot(rij, dv1j_phi) / r * Ajy, -1.*np.dot(rij, dv1j_phi) / r * Ajz]
#             # Just Missing a kb term now
#
#             #Use Relationships to calc Hii, Hjj
#             PA_ii = -1*PA_ji
#             AP_ii = -1*AP_ij
#
#             PA_jj = -1*PA_ij
#             AP_jj = -1*AP_ji
#
#             # Hii Angle to Angle bending terms
#             AAB_ii = np.full((2, 2), 0.0)
#             AAB_ii[0] = [np.dot(rij, dv1i_theta)**2. / (r**2.), np.dot(rij, dv1i_theta) * np.dot(rij, dv1i_phi) / (r**2.)]
#             AAB_ii[1] = [np.dot(rij, dv1i_theta) * np.dot(rij, dv1i_phi) / (r**2.), np.dot(rij, dv1i_phi)**2. / (r**2.)]
#             # Just missing kb term
#
#             # Hjj Angle to Angle bending Terms
#             AAB_jj = np.full((2, 2), 0.0)
#             AAB_jj[0] = [np.dot(rij, dv1j_theta)**2. / (r**2.), np.dot(rij, dv1j_theta) * np.dot(rij, dv1j_phi) / (r**2.)]
#             AAB_jj[1] = [np.dot(rij, dv1j_theta) * np.dot(rij, dv1j_phi) / (r**2.), np.dot(rij, dv1j_phi)**2. / (r**2.)]
#             # Just missing kb term
#
#             #Hii Angle to Angle torsion Terms
#             AAT_ii = np.full((2, 2), 0.0)
#             AAT_ii[0] = [np.dot(b1, dv1i_theta)**2., np.dot(b1, dv1i_phi) * np.dot(b1, dv1i_theta)]
#             AAT_ii[1] = [np.dot(b1, dv1i_phi) * np.dot(b1, dv1i_theta), np.dot(b1, dv1i_phi)**2.]
#
#             #Hjj Angle to Angle torsion Terms
#             AAT_jj = np.full((2, 2), 0.0)
#             AAT_jj[0] = [np.dot(a1, dv1j_theta)**2., np.dot(a1, dv1j_phi) * np.dot(a1, dv1j_theta)]
#             AAT_jj[1] = [np.dot(a1, dv1j_phi) * np.dot(a1, dv1j_theta), np.dot(a1, dv1j_phi)**2.]
#
#             # Same for ii, jj, and ij terms
#             # Torsion Derviatives a1 & b1
#             AAT_ij = np.full((2, 2), 0.0)
#             AAT_ij[0] = [np.dot(a1, dv1j_theta) * np.dot(b1, dv1i_theta), np.dot(a1, dv1j_phi) * np.dot(b1, dv1i_theta)]
#             AAT_ij[1] = [np.dot(a1, dv1j_theta) * np.dot(b1, dv1i_phi), np.dot(a1, dv1j_phi) * np.dot(b1, dv1i_phi)]
#             # Just missing kt term now
#
#             # Torsion Derivatives ai, b1 for Hji
#             AAT_ji = np.matrix.transpose(AAT_ij)
#             # AAT_ji = AAT_ij
#
#             # Torsion Derviatives a3 & b3
#             AA3T_ij = np.full((2, 2), 0.0)
#             AA3T_ij[0] = [np.dot(a3, dv3j_theta) * np.dot(b3, dv3i_theta), np.dot(a3, dv3j_phi) * np.dot(b3, dv3i_theta)]
#             AA3T_ij[1] = [np.dot(a3, dv3j_theta) * np.dot(b3, dv3i_phi), np.dot(a3, dv3j_phi) * np.dot(b3, dv3i_phi)]
#             # just missing kt term now
#
#             # Torsion Derivatives a3, b3 for Hji
#             AA3T_ji = np.matrix.transpose(AA3T_ij)
#             # AA3T_ji = AA3T_ij
#
#             # Hii Angle to Angle torsion Terms a3 & b3
#             AA3T_ii = np.full((2, 2), 0.0)
#             AA3T_ii[0] = [np.dot(b3, dv3i_theta) ** 2., np.dot(b3, dv3i_phi) * np.dot(b3, dv3i_theta)]
#             AA3T_ii[1] = [np.dot(b3, dv3i_phi) * np.dot(b3, dv3i_theta), np.dot(b3, dv3i_phi) ** 2.]
#
#             # Hjj Angle to Angle torsion Terms a3& b3
#             AA3T_jj = np.full((2, 2), 0.0)
#             AA3T_jj[0] = [np.dot(a3, dv3j_theta) ** 2., np.dot(a3, dv3j_phi) * np.dot(a3, dv3j_theta)]
#             AA3T_jj[1] = [np.dot(a3, dv3j_phi) * np.dot(a3, dv3j_theta), np.dot(a3, dv3j_phi) ** 2.]
#
#             #All derivatives Calculated, Now just Fill Hij
#             subHij = np.full((7, 7), 0.0)
#             subHij[0:3, 0:3] = -dir * kb
#             subHij[0:3, 3:5] = kb * PA_ij
#             subHij[3:5, 0:3] = kb * AP_ij
#             subHij[3:5, 3:5] = kt * AAT_ij
#             subHij[5:7, 5:7] = kt * AA3T_ij
#
#             #Hji
#             subHji = np.full((7, 7), 0.0)
#             subHji[0:3, 0:3] = -dir * kb
#             subHji[0:3, 3:5] = kb * PA_ji
#             subHji[3:5, 0:3] = kb * AP_ji
#             subHji[3:5, 3:5] = kt * AAT_ji
#             subHji[5:7, 5:7] = kt * AA3T_ji
#
#             #Fill Hii
#             subHii = np.full((7, 7), 0.0)
#             subHii[0:3, 0:3] = dir * kb
#             subHii[0:3, 3:5] = kb * PA_ii
#             subHii[3:5, 0:3] = kb * AP_ii
#             subHii[3:5, 3:5] = kt * AAT_ii + kb * AAB_ii
#             subHii[5:7, 5:7] = kt * AA3T_ii
#
#             #Fill Hjj
#             subHjj = np.full((7, 7), 0.0)
#             subHjj[0:3, 0:3] = dir * kb
#             subHjj[0:3, 3:5] = kb * PA_jj
#             subHjj[3:5, 0:3] = kb * AP_jj
#             subHjj[3:5, 3:5] = kt * AAT_jj + kb * AAB_jj
#             subHjj[5:7, 5:7] = kt * AA3T_jj
#
#             #Fill Actual Matrix
#             bend_hess[7 * i:7 * i + 7, 7 * j:7 * j + 7] += subHij
#             bend_hess[7 * j:7 * j + 7, 7 * i:7 * i + 7] += subHji
#             bend_hess[7 * i:7 * i + 7, 7 * i:7 * i + 7] += subHii
#             bend_hess[7 * j:7 * j + 7, 7 * j:7 * j + 7] += subHjj
#
#         return bend_hess
#
#     def total_hess(self, bend_hess, spring_hess, ks=1):
#         for i in range(self.cc):
#             for j in range(self.cc):
#                 bend_hess[7*i:7*i+3, 7*j:7*j+3] += spring_hess[3*i:3*i+3, 3*j:3*j+3]*ks
#         return bend_hess
#
#     def calc_bfacts(self, inversehess):
#         msd = []
#         for i in range(self.cc):
#             print(inversehess[7*i:7*i+3, 7*i:7*i+3])
#             msd.append(self.kb*self.T*np.sum(np.diag(inversehess[7*i:7*i+3, 7*i:7*i+3])))
#             #this ain't it
#             # msd.append(self.kb*self.T*np.sum(np.diag(inversehess[7*i:7*i+7, 7*i:7*i+7])))
#         ana_bfactors = [self.bconv * x for x in msd]
#
#         return ana_bfactors


#Known Issues: Will not converge if under-constrained residues present, will have to raise cutoff value
class HANM(ANM):
    def __init__(self, coord, exp_bfactors, cutoff=15, T=300, scale_factor=0.3, mcycles=5, ncycles=7):
        super().__init__(coord, exp_bfactors, T=T, cutoff=cutoff)
        self.spring_constant_matrix = np.full((self.cc, self.cc), 1.)
        self.calc_ANM_unitary(cuda=False)
        self.spring_constant_matrix = np.full((self.cc, self.cc), self.ana_gamma)

        self.restraint_force_constants = []
        self.bond_fluctuations = []
        self.bond_fluctuations0 = []

        self.scale_factor = scale_factor
        self.mcycles = mcycles
        self.ncycles = ncycles
        self.routine_finished = False
        self.model_id = 'HANM'

    def hanm_calc_restraint_force_constant(self, bcal):
        restraint_force_constants = []
        for i in range(self.cc):
            ki_res = self.scale_factor * self.kb * self.T * 8 * math.pi ** 2. * (bcal[i] - self.exp_bfactors[i]) / \
                     (bcal[i] * self.exp_bfactors[i])
            restraint_force_constants.append(ki_res)
        # print(restraint_force_constants)
        return restraint_force_constants

    def hanm_add_restraints(self, hess, restraint_force_constants):
        for i in range(self.cc):
            hess[3 * i, 3 * i] += restraint_force_constants[i]
            hess[3 * i + 1, 3 * i + 1] += restraint_force_constants[i]
            hess[3 * i + 2, 3 * i + 2] += restraint_force_constants[i]

    def hanm_calc_bond_fluctuations(self, hess, cuda=False):
        if cuda:
            thess = np.asarray(hess, dtype=np.float32, order='C')
            cu_hess = gpuarray.to_gpu(thess)
            cu_evecs, cu_evals = cuda_la.eig(cu_hess, imag='T')
            evals, evecs = cu_evals.get(), cu_evecs.get()
        else:
            evals, evecs = la.eig(hess)
        # print(self.evals)
        idx = np.argsort(abs(evals))
        evals = np.asarray(evals[idx])

        if cuda:
            evecs = np.asarray(evecs[idx, :])
            evecs = np.swapaxes(evecs, 1, 0)
        else:
            evecs = np.asarray(evecs[:, idx])

        # fig = plt.figure()
        # plt.imshow(evecs)
        # plt.savefig('evecs.png', dpi=600)

        # print('eVALS:', evals[6], evals[7])
        # print('evecs:', evecs[0, 6], evecs[0, 7], evecs[1, 6], evecs[1, 7])

        bcal = [0. for x in range(self.cc)]
        bond_fluc = np.full((self.cc, self.cc), 0.)
        for i in range(self.cc):
            for j in range(6, 3 * self.cc):
                if evals[j] != 0.:
                    bcal[i] += np.inner(evecs[3 * i: 3 * i + 3, j], evecs[3 * i: 3 * i + 3, j]) / evals[j]
            bcal[i] *= self.kb * self.T * self.bconv
        for i in range(self.cc):
            for j in range(self.cc):
                dis = self.distance_matrix[i, 4 * j]
                if dis:
                    tmp = np.asarray(self.distance_matrix[i, 4 * j + 1:4 * j + 4])
                    # print(tmp)
                    for k in range(6, 3 * self.cc):
                        p = tmp / dis * (evecs[3 * i:3 * i + 3, k] - evecs[3 * j:3 * j + 3, k])
                        if evals[k] != 0.:
                            bond_fluc[i, j] += np.sum(p) ** 2. / evals[k]
        return bcal, bond_fluc

    def hanm_nma(self, fc, fc_res, cuda=False):
        hess = self.calc_hess_fast_sc(fc)
        self.hanm_add_restraints(hess, fc_res)
        bcal, bond_fluc = self.hanm_calc_bond_fluctuations(hess, cuda=cuda)
        return bcal, bond_fluc

    def routine(self, cuda=False):

        mthreshold1 = 0.005  # relative
        mthreshold2 = 0.01  # absolute
        nthreshold = 0.001
        alpha = 1.0
        bcal = []
        bcalprev = [0. for x in range(self.cc)]

        fc_res0 = [0. for x in range(self.cc)]
        sc_mat_tmp = np.full((self.cc, self.cc), self.ana_gamma)

        bcal, bond_fluc = self.hanm_nma(sc_mat_tmp, fc_res0, cuda=cuda)

        for i in range(self.mcycles):

            mcheck = 1
            for y in range(self.cc):
                rb1 = abs(bcal[y] - bcalprev[y]) / bcal[y]
                rb2 = abs(bcal[y] - self.exp_bfactors[y]) / self.exp_bfactors[y]

                # criterion can be changed if need be
                if rb2 > mthreshold2:
                    mcheck = 0.

            bcalprev = [x for x in bcal]

            if mcheck:
                self.ana_bfactors = bcal
                self.spring_constant_matrix = sc_mat_tmp
                self.routine_finished = True
                print("Bfactors have converged")
                break

            # Calc Restraint forces, add to hessian, calc bond_fluc, save as bond_fluctuations0
            fc_res = self.hanm_calc_restraint_force_constant(bcal)

            bcal, bond_fluc0 = self.hanm_nma(sc_mat_tmp, fc_res, cuda=cuda)


            for j in range(self.ncycles):

                bcal, bond_fluc = self.hanm_nma(sc_mat_tmp, fc_res0, cuda=cuda)

                ncheck = 1
                for x in range(self.cc):
                    for y in range(self.cc):
                        if x >= y:
                            continue
                        if self.distance_matrix[x, 4 * y] != 0.:
                            delta_fluc = bond_fluc[x, y] - bond_fluc0[x, y]
                            r2 = abs(delta_fluc) / bond_fluc0[x, y]
                            if r2 > nthreshold:
                                ncheck = 0.

                        # update Spring Constant Matrix
                        sc_mat_tmp[x, y] = 1. / ((1. / sc_mat_tmp[x, y]) - alpha * delta_fluc)
                        sc_mat_tmp[y, x] = sc_mat_tmp[x, y]

                if ncheck:
                    print("Force Constants have converged after %d mcycles and %d ncycles" % (i, j))
                    self.ana_bfactors = bcal
                    self.spring_constant_matrix = sc_mat_tmp
                    self.routine_finished = True
                    break
                else:
                    print("Force Constants have not converged after %d mcycles and %d ncycles" % (i, j))
        self.ana_bfactors = bcal
        self.spring_constant_matrix = sc_mat_tmp
        self.routine_finished = True

    def hanm_theor_bfactors(self, outfile):
        if self.routine_finished:
            free_compare(outfile, self.exp_bfactors, self.ana_bfactors,
                         legends=['Experimental  (PDB)', 'Analytical (HANM)'])
        else:
            print('HANM has not been run')


#Helper function for getting coordinates
def find_nearest(array,value):
    array=np.asarray(array)
    idx=(np.abs(array - value)).argmin()
    return array[idx]
#Easier to write
s=' '
n='\n'


#Giant File Writer
class protein:
    def __init__(self, pdbfile, cutoff=15, pottype='s', potential=5.0, offset_indx=0, strand_offset=0, backbone_weight=0,
                 importscmatrix=False, scmatrix=0, angconstrain=False, upstreamdir=''):
        self.su = 8.518
        self.boxsize = 0
        self.pi = 0
        self.backbone_weight = backbone_weight
        self.angconstrain = angconstrain

        if "/" in pdbfile:
            pdbid = pdbfile.rsplit('/', 1)[1].split('.')[0]
        else:
            pdbid = pdbfile.split('.')[0]

        wdir = os.getcwd()
        # Outfiles they go into directory you call script from
        self.parfile = wdir + os.path.join(upstreamdir + '/generated.par')
        self.topfile = wdir + os.path.join(upstreamdir + '/generated.top')
        self.datfile = wdir + os.path.join(upstreamdir + '/generated.dat')

        self.sim_force_const = .05709
        self.pottype = pottype
        self.import_sc = importscmatrix

        if importscmatrix:
            self.spring_constant_matrix = scmatrix
            self.potential = 0.
        else:
            self.spring_constant_matrix = []
            self.potential = potential

        self.rc = cutoff
        self.topbonds = []
        self.chainnum = 0

        self.strand_offset = strand_offset
        self.offset_indx = offset_indx
        self.topology = []
        self.conv = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P', 'THR': 'T',
                     'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'ALA': 'A',
                     'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

        self.a1s = []
        self.a3s = []
        self.getalphacarbons(pdbid, pdbfile)

    def three_to_one_seq(self, res):
        ol = self.conv[res]
        return ol

    def adjust_coord(self, column, min):
        for i in range(0, self.pi):
            if min > 0:
                self.coord[i, column] -= min
            else:
                self.coord[i, column] += min

    def getalphacarbons(self, pdbid, pdbfile):
        structure = Bio.PDB.PDBParser().get_structure(pdbid, pdbfile)
        model = Bio.PDB.Selection.unfold_entities(structure, 'C')

        # orientation vector container
        tmp_N_vectors = []
        # chain index, residue index, residue identitiy, CA coordinates
        cindx, rindx, rids, coord = [], [], [], []
        if len(model) == 1:
            self.single_chain = True
        else:
            self.single_chain = False
        for cid, chain in enumerate(model):
            for rid, residue in enumerate(chain.get_residues()):
                tags = residue.get_full_id()
                if tags[3][0] == " ":
                    onelettercode = self.three_to_one_seq(residue.get_resname())
                    atoms = residue.get_atoms()
                    # Center of Mass Used for orientations
                    com = np.full(3, 0.0)
                    count = 0
                    for atom in atoms:
                        if atom.get_id() == 'CA':
                            coordinates = atom.get_coord()
                            cindx.append(cid)
                            rindx.append(rid)
                            rids.append(onelettercode)
                            coord.append(coordinates)
                        if atom.get_id() != 'N' or atom.get_id() != 'C' or atom.get_id() != 'O':
                            count += 1
                            com += atom.get_coord()
                    if np.sum(com) != 0.:
                        com /= count
                        # Vector from CA coordinates to Center of Mass
                        nvec = normalize_vector(com - coord[-1])
                        tmp_N_vectors.append(nvec)
                    else:
                        tmp_N_vectors.append(np.asarray([0., 0., 0.]))
        self.topology = zip(cindx, rindx, rids)
        acs = np.divide(np.array(coord), self.su)
        self.coord = acs

        # Dealing with Null Normal Vectors
        normal_vectors = []
        fcoord = np.asarray(flatten(self.coord))
        for iid, i in enumerate(tmp_N_vectors):
            if np.sum(i ** 2) == 0.:
                # vector from particle i to i + 1
                rij = fcoord[iid + 1] - fcoord[iid]
                # vector from particle i to i-1
                rih = fcoord[iid - 1] - fcoord[iid]
                nvec = normalize_vector(np.cross(rij, rih))
                normal_vectors.append(nvec)
            else:
                nvec = normalize_vector(i)
                normal_vectors.append(nvec)
        normal_vectors = np.asarray(normal_vectors)

        #Calculated a1 vectors
        self.a1s = normal_vectors

        #Calc a3 vectors
        basis_vecs = np.asarray([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
        # Threshold to declare whether vector is perpendicular or parallel
        epsilon = 1e-3
        tmp_orth_vectors = []
        for a1 in self.a1s:
            if np.linalg.norm(a1) != 1.:
                a1 = normalize_vector(a1)
            # check if parallel to 1st basis vector
            if not (np.dot(a1, basis_vecs[0]) > 1 - epsilon):
                a3 = np.cross(a1, basis_vecs[0])
            elif not (np.dot(a1, basis_vecs[1]) > 1 - epsilon):
                a3 = np.cross(a1, basis_vecs[1])
            # Arbitrarily Defined Just needs to be perpendicular
            tmp_orth_vectors.append(a3)
        self.a3s = np.asarray(tmp_orth_vectors)


    def WriteParFile_custombonds(self, bondlist, potential=5.0):
        make = open(self.parfile, 'w')
        self.pi = len(self.coord)
        make.write(str(len(self.coord)))
        make.write(n)
        make.close()
        p = open(self.parfile, "a")
        custom_potentials = len(bondlist[0]) > 2
        ordered = sorted(bondlist, key=lambda tup: (tup[0], tup[1]))

        if custom_potentials:
            for bond in ordered:
                i, j, po = bond
                if j-i == 1:
                    dx = self.coord[j, 0] - self.coord[i, 0]
                    dy = self.coord[j, 1] - self.coord[i, 1]
                    dz = self.coord[j, 2] - self.coord[i, 2]
                    dist = math.sqrt(abs(dx) ** 2 + abs(dy) ** 2 + abs(dz) ** 2)
                    r = np.asarray([dx / dist, dy / dist, dz / dist])
                    a0, b0, c0, d0 = np.dot(self.a1s[i], r), np.dot(self.a1s[j], -1. * r), np.dot(self.a1s[i], self.a1s[j]), np.dot(self.a3s[i], self.a3s[j])
                    print(i + self.offset_indx, j + self.offset_indx, dist, 's', po+self.backbone_weight,
                          a0, b0, c0, d0,
                          file=p)
                else:
                    self.topbonds.append((i, j))
                    dx = self.coord[j, 0] - self.coord[i, 0]
                    dy = self.coord[j, 1] - self.coord[i, 1]
                    dz = self.coord[j, 2] - self.coord[i, 2]
                    dist = math.sqrt(abs(dx) ** 2 + abs(dy) ** 2 + abs(dz) ** 2)
                    print(i + self.offset_indx, j + self.offset_indx, dist, 's', p, file=p)
        else:
            for bond in ordered:
                i, j = bond
                if j - i == 1:
                    dx = self.coord[j, 0] - self.coord[i, 0]
                    dy = self.coord[j, 1] - self.coord[i, 1]
                    dz = self.coord[j, 2] - self.coord[i, 2]
                    dist = math.sqrt(abs(dx) ** 2 + abs(dy) ** 2 + abs(dz) ** 2)
                    r = np.asarray([dx / dist, dy / dist, dz / dist])
                    a0, b0, c0, d0 = np.dot(self.a1s[i], r), np.dot(self.a1s[j], -1. * r), np.dot(self.a1s[i],
                                                                                                  self.a1s[j]), np.dot(
                        self.a3s[i], self.a3s[j])
                    print(i + self.offset_indx, j + self.offset_indx, dist, 's', potential + self.backbone_weight,
                          a0, b0, c0, d0,
                          file=p)
                else:
                    self.topbonds.append((i, j))
                    dx = self.coord[j, 0] - self.coord[i, 0]
                    dy = self.coord[j, 1] - self.coord[i, 1]
                    dz = self.coord[j, 2] - self.coord[i, 2]
                    dist = math.sqrt(abs(dx) ** 2 + abs(dy) ** 2 + abs(dz) ** 2)
                    print(i + self.offset_indx, j + self.offset_indx, dist, 's', potential, file=p)

        p.close()


    def WriteParFile(self, pottype='s', potential=5.0):
        make = open(self.parfile, 'w')
        self.pi = len(self.coord)  # Particle Index for range operations
        #print(self.coord.shape, len(self.coord))
        print('N = ', self.pi)
        make.write(str(len(self.coord)))
        make.write(n)
        make.close()
        p = open(self.parfile, "a")
        for i in range(0, self.pi):
            for j in range(0, self.pi):
                if i >= j:
                    continue
                else:
                    dx = self.coord[j, 0] - self.coord[i, 0]
                    dy = self.coord[j, 1] - self.coord[i, 1]
                    dz = self.coord[j, 2] - self.coord[i, 2]
                    dist = math.sqrt(abs(dx) ** 2 + abs(dy) ** 2 + abs(dz) ** 2)
                    if self.pottype == 's':
                        if dist < (self.rc / self.su):
                            self.topbonds.append((i, j))
                            if abs(i - j) == 1:
                                if self.import_sc:
                                    if self.angconstrain:
                                        r = np.asarray([dx/dist, dy/dist, dz/dist])
                                        a0, b0, c0, d0 = np.dot(self.a1s[i], r), np.dot(self.a1s[j], -1.*r), np.dot(self.a1s[i], self.a1s[j]), np.dot(self.a3s[i], self.a3s[j])
                                        spring_constant = self.spring_constant_matrix[i, j] / self.sim_force_const
                                        print(i + self.offset_indx, j + self.offset_indx, dist, self.pottype, spring_constant,
                                              a0, b0, c0, d0,
                                              file=p)
                                    else:
                                        spring_constant = self.spring_constant_matrix[i, j] / self.sim_force_const
                                        print(i + self.offset_indx, j + self.offset_indx, dist, self.pottype, spring_constant,
                                              file=p)
                                else:
                                    if self.angconstrain:
                                        r = np.asarray([dx/dist, dy/dist, dz/dist])
                                        a0, b0, c0, d0 = np.dot(self.a1s[i], r), np.dot(self.a1s[j], -1.*r), np.dot(self.a1s[i], self.a1s[j]), np.dot(self.a3s[i], self.a3s[j])
                                        print(i + self.offset_indx, j + self.offset_indx, dist, self.pottype, self.potential,
                                              a0, b0, c0, d0,
                                              file=p)
                                    else:
                                        print(i + self.offset_indx, j + self.offset_indx, dist, self.pottype,
                                              self.potential, file=p)
                            else:
                                if self.import_sc:
                                    spring_constant = self.spring_constant_matrix[i, j] / self.sim_force_const
                                    print(i + self.offset_indx, j + self.offset_indx, dist, self.pottype, spring_constant,
                                          file=p)
                                else:
                                    print(i + self.offset_indx, j + self.offset_indx, dist, self.pottype, self.potential,
                                          file=p)
        p.close()

    def WriteConfFile(self):
        xcoord = np.array(self.coord[:-1, 0])
        ycoord = np.array(self.coord[:-1, 1])
        zcoord = np.array(self.coord[:-1, 2])
        xmin = find_nearest(xcoord, 0)
        ymin = find_nearest(ycoord, 0)
        zmin = find_nearest(zcoord, 0)
        self.adjust_coord(0, xmin)
        self.adjust_coord(1, ymin)
        self.adjust_coord(2, zmin)
        span = np.array(
            [np.max(xcoord) - np.min(xcoord), np.max(ycoord) - np.min(ycoord), np.max(zcoord) - np.min(zcoord)])
        self.boxsize = 2.5 * math.ceil(np.max(span))
        conf = open(self.datfile, "w")
        print('t = 0', file=conf)
        print("b =", str(self.boxsize), str(self.boxsize), str(self.boxsize), file=conf)
        print("E =", 0, 0, 0, file=conf)
        for i in range(0, self.pi):
            print(self.coord[i, 0], self.coord[i, 1], self.coord[i, 2], ' '.join(map(str, self.a1s[i])), ' '.join(map(str, self.a3s[i])), 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, file=conf)
        conf.close()

    def WriteTopFile(self):
        cindx, rindx, rids = zip(*self.topology)
        # Finding where the chains start and stop
        t_aft = [cindx.index(x) for x in np.arange(1, cindx[-1] + 1, 1)]
        t_pri = [x - 1 for x in t_aft]

        t = open(self.topfile, 'w')
        print(self.pi, len(t_aft) + 1, file=t)
        if self.pottype == 's':
            # Get bonds
            fullneighs = []
            for j in range(0, self.pi):
                neighs = []
                for x, y in self.topbonds:
                    if x + 1 != y and x == j:
                        neighs.append(y)
                nebs = list(neighs)
                nebs_adj = [x + self.offset_indx for x in nebs]
                fullneighs.append(nebs_adj)

            for cid, i in enumerate(cindx):
                ci_adj = -1 * (i + self.strand_offset + 1)
                rindx_adj = cid + self.offset_indx
                olc = rids[cid]
                bonds = fullneighs[cid]
                #print(bonds)
                if self.pottype == 's':
                    if cid == self.pi - 1 or cid in t_pri:
                        print(ci_adj, olc, rindx_adj - 1, -1, *bonds, file=t)
                    elif cid == 0 or cid in t_aft:
                        print(ci_adj, olc, -1, rindx_adj + 1, *bonds, file=t)
                    else:
                        print(ci_adj, olc, rindx_adj - 1, rindx_adj + 1, *bonds, file=t)
        t.close()

    def WriteSimFiles(self):
        self.WriteParFile()
        print('Wrote Par File to', self.parfile)
        self.WriteConfFile()
        print('Wrote Configuration (dat) File to', self.datfile)
        self.WriteTopFile()
        print('Wrote Topology File to', self.topfile)



def export_to_simulation(model, pdbfile, upstreamdir = ''):
    if model.model_id == 'ANM':
        p = protein(pdbfile, cutoff=model.cutoff, potential=model.ana_gamma/0.05709, upstreamdir=upstreamdir)
        p.WriteSimFiles()
    elif model.model_id == 'HANM' or model.model_id == 'MVP' or model.model_id == 'mANM':
        p = protein(pdbfile, cutoff=model.cutoff, importscmatrix=True, scmatrix=model.spring_constant_matrix, upstreamdir=upstreamdir)
        p.WriteSimFiles()
    elif model.model_id == 'pep':
        p = protein(pdbfile, cutoff=0, potential=0, angconstrain=True, backbone_weight=model.backbone_weight, upstreamdir=upstreamdir)
        p.WriteParFile_custombonds(model.bonds, potential=model.potential)
        p.WriteConfFile()
        p.WriteTopFile()
    elif model.model_id == 'ANMT':
        p = protein(pdbfile, cutoff=model.cutoff, potential=model.ana_gamma/0.05709, angconstrain=True, upstreamdir=upstreamdir)
        p.WriteSimFiles()
    else:
        print('Model Type', model.model_id, 'is not supported')
