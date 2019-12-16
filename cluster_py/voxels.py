#!/usr/bin/env python
# encoding: utf-8
# author : Lianjun
# created time: 2019-12-12

from math import sqrt, exp
import numpy as np
import pymp

from rdkit import Chem
import moleculekit.smallmol.smallmol as mk
import moleculekit.tools.voxeldescriptors as mk_voxel
from moleculekit.util import uniformRandomRotation

NVOXEL_X = 32  # can be tuned
NVOXEL_Y = NVOXEL_X  # can be tuned
NVOXEL_Z = NVOXEL_X  # can be tuned
VOXELSIZE = 1.0  # can be tuned
# N_VOXEL's are the number of grids per axis.
# VOXELSIZE, i.e., the resolution, has unit Angstrom.
# NVOXEL_X * VOXELSIZE is the length (A) of the box at the x-axis.
# Make sure the box is big enough to hold the molecule.
CHANNELS = {"Hydrophobic": 0,
            "Aromatic": 1,
            "Acceptor": 2,
            "Donor": 3,
            "HeavyAtom": 4,
            "Pharm_Aromatic": 5,
            "Pharm_Acceptor": 6,
            "Pharm_Donor": 7}
ADD_PHARM = True  # can be tuned
# if ADD_PHARM = False, then the following channels are used
CHANNELS_INIT = {"Hydrophobic": 0,
                 "Aromatic": 1,
                 "Acceptor": 2,
                 "Donor": 3,
                 "PosIonizable": 4,
                 "NegIonizable": 5,
                 "Metal": 6,
                 "HeavyAtom": 7}
NCHANNELS = 8
NCHANNELS_MOL = 5
NCHANNELS_PHARM = 3

# voxelization
PHARM_RADIUS = 2.0
DISPLACEMENT = 4
IS_SHUFFLE = False
TRUNC = 5
N_PROC = 4  # number of processors


class VoxelMol:
    """
    manipulate molecules with conformations and voxel representations

    Attributes
    -------------
    mol : rdkit.Chem.rdchem.Mol object
    voxels : 4D numpy array [NVOXEL_X][NVOXEL_Y][NVOXEL_Z][NCHANNELS]

    Parameters
    -------------
    inputmol : str or rdkit.Chem.rdchem.Mol object
        Path to mol2 file, or Mol object

    Examples
    -------------
    >>> from voxels import VoxelMol

    # create a VoxelMol object from mol2 file
    >>> vmol = VoxelMol("example.mol2")
    # or create a VoxelMol object from rdkit.Mol
    >>> rdkitmol = vmol.mol
    >>> vmol = VoxelMol(rdkitmol)

    # use "vmol.voxels" to obtain a numpy 4D array (x, y, z, channel)
    >>> vmol.voxels.shape
    (32, 32, 32, 8)

    # one may flatten it into 1D array
    >>> flattened = vmol.voxels.flatten()
    array([ 0.,  0.,  0., ...,  0.,  0.,  0.], dtype=float32)
    >>> flattened.shape
    (262144,)
    # or select channels of interest
    >>> vmol.voxels[:,:,:,[0,4]].shape
    (65536,)

    # if VMD is installed, one may visualize the voxel points by channel
    >>> vmol.view_voxels()

    # see how many voxel points have non-zero values at each channel
    >>> vmol.num_grids_per_channel()
    {'Hydrophobic': 877,
     'Aromatic': 877,
     'Acceptor': 0,
     'Donor': 0,
     'HeavyAtom': 877,
     'Pharm_Aromatic': 494,
     'Pharm_Acceptor': 0,
     'Pharm_Donor': 0}

    # obtain voxel points ([x, y, z]) that have value at the given channel
    >>> vmol.get_index_by_channel(1) # or vmol.get_index_by_channel("Aromatic")
    array([[10, 14, 16],
           [10, 15, 15],
           [10, 15, 16],
           ...,
           [22, 17, 16],
           [22, 17, 17],
           [22, 18, 16]])
    """

    def __init__(self, inputmol):
        if isinstance(inputmol, str):
            self.mol = Chem.MolFromMol2File(inputmol)
        elif isinstance(inputmol, Chem.rdchem.Mol):
            self.mol = inputmol
        else:
            raise AttributeError("Please input filepath-to-mol2 or rdkit mol")
        self.voxels = self.voxelize()

        # if isinstance(inputmol, str):
        #     self.mol = Chem.SDMolSupplier(inputmol)
        # elif isinstance(inputmol, Chem.rdchem.Mol):
        #     self.mol = inputmol
        # else:
        #     raise AttributeError("Please input filepath-to-mol2 or rdkit mol")
        # self.voxels = self.voxelize()

    def _shuffle_loc(self, coords, by_range=1.0):
        """
        This function randomly re-locates the coordinates within a range
        Only called when IS_SHUFFLE = True

        Parameters
        -------------
        coords : numpy array
            for molecule coordinates: 2D numpy array of shape(#items, 3)
            for molecule center: 1D numpy array [x, y, z]
        by_range : float, default +/- 0.5

        Returns
        -------------
        coords_shuffled : numpy array

        Examples
        -------------
        >>> voxmol._shuffle_loc(np.array([[0,0,0],[1,1,1]]))
        array([[ 0.3391448 ,  0.30718527, -0.00914312],
               [ 1.19678919,  0.54693674,  1.41423   ]])
        >>> voxmol._shuffle_loc(np.array([0,0,0]), 4)
        array([-1.92925533,  1.93960189,  1.65201393])
        """
        return coords + (np.random.rand(*coords.shape) - 0.5) * by_range

    def _get_aromatic_groups(self, rdkit_mol):
        """
        This function returns a list of aromatic rings containing atom indices

        Parameters
        -------------
        rdkit_mol : rdkit.Chem.rdchem.Mol object with one conformation

        Returns
        -------------
        groups : a list of tuple
            each tuple is a group of atom indices that belong to one
            aromatic ring

        Examples
        -------------
        >>> voxmol = VoxelMol('c1ccccc1')
        >>> voxmol._get_aromatic_groups(voxmol.mol)
        [(0, 5, 4, 3, 2, 1)]
        """

        groups = []
        ring_atoms_tuple = rdkit_mol.GetRingInfo().AtomRings()
        for ring_group in ring_atoms_tuple:
            if all([rdkit_mol.GetAtomWithIdx(x).GetIsAromatic() for x in
                    ring_group]):
                groups.append(ring_group)
        return groups

    def _get_aromatic_centers(self, coords, aromatic_groups):
        """
        This function returns the coordinates of the center of each ring

        Parameters
        -------------
        coords : 2D numpy array of shape(#atoms, 3)
        aromatic_groups : a list of tuple
            each tuple is a group of atom indices that belong to one
            aromatic ring

        Returns
        -------------
        aromatic_centers : 2D numpy array of shape(#rings, 3)

        Examples
        -------------
        >>> voxmol = VoxelMol('c1ccccc1')
        >>> voxmol._get_aromatic_centers(
            voxmol.get_coords(), voxmol._get_aromatic_groups(voxmol.mol))
        array([[ -1.09883584e-08,  -1.00937274e-07,  -9.66406850e-07]])
        """

        if len(coords) == 0 or len(aromatic_groups) == 0:
            return np.empty(shape=(0, 3))
        aromatic_centers = np.array([coords[np.array(a_group)].mean(axis=0) for
                                     a_group in aromatic_groups])
        if len(aromatic_centers) == 0:  # Make sure the shape is correct
            aromatic_centers = aromatic_centers.reshape(
                aromatic_centers.shape[0], 3)
        return aromatic_centers

    def add_pharmacophore(self, mol_coords, channel_init, rdkit_mol):
        """
        This function appends the coordinates and channels of the
        pharmocophores of the molecule to the molecular ones

        Parameters
        -------------
        mol_coords : 2D numpy array of shape(#atoms, 3)
        channel_init : 2D numpy array of shape(#atoms, 8)
        rdkit_mol : rdkit.Chem.rdchem.Mol object with one conformation

        Returns
        -------------
        coords_all : 2D numpy array of shape((#atoms + #aromatic + #acceptor +
                                             #donor), 3)
        channel_all : 2D numpy array of shape((#atoms + #aromatic + #acceptor +
                                             #donor), 8)

        Examples
        -------------
        >>> voxmol = VoxelMol('CNC(=O)CCCN=c1[nH]c(=N)[nH]c(=N)[nH]1')
        # inside voxelize()
        >>> print(self.add_pharmacophore(coords_mol, channel_mol, rdkit_mol))
        (array([[-5.28952219, -0.37745559,  1.21308548],
                [-4.47566669, -0.73377859,  0.07848694],
                [-3.29904781, -0.06621711, -0.19182064],
                [-2.90531088,  0.89467319,  0.45775558],
                [-2.55246381, -0.65415994, -1.37476359],
                [-1.24934298,  0.06868279, -1.70766534],
                [-0.16139219, -0.14784289, -0.6543105 ],
                [ 1.09500838,  0.38192878, -1.18516229],
                [ 2.11657184,  0.40037646, -0.39792888],
                [ 3.31152969,  0.89155036, -0.86571289],
                [ 4.4958691 ,  0.98840379, -0.19373445],
                [ 5.59747835,  1.45881372, -0.67236686],
                [ 4.41886294,  0.52369193,  1.08816941],
                [ 3.31993959,  0.00723685,  1.71181562],
                [ 3.2784688 , -0.42392485,  2.92681836],
                [ 2.21702067, -0.01915432,  0.90657477],
                [ 3.59555839,  0.91522893,  0.8531397 ],
                [-2.93992768,  0.43075294,  0.11271633],
                [-4.42957192, -1.00632129,  0.54498008],
                [ 3.44299446,  1.34659459, -0.77300089],
                [ 5.82812742,  1.227013  , -0.59920996],
                [ 4.91605467,  0.48422669,  1.36305552],
                [ 3.24463503, -0.15482432,  2.42722583],
                [ 2.33468746, -0.26546558,  1.11435183]]),
         array([[ 0.  ,  0.  ,  0.  ,  0.  ,  1.7 ,  0.  ,  0.  ,  0.  ],
                [ 0.  ,  0.  ,  0.  ,  1.55,  1.55,  0.  ,  0.  ,  0.  ],
                [ 0.  ,  0.  ,  0.  ,  0.  ,  1.7 ,  0.  ,  0.  ,  0.  ],
                [ 0.  ,  0.  ,  1.52,  0.  ,  1.52,  0.  ,  0.  ,  0.  ],
                [ 1.7 ,  0.  ,  0.  ,  0.  ,  1.7 ,  0.  ,  0.  ,  0.  ],
                [ 1.7 ,  0.  ,  0.  ,  0.  ,  1.7 ,  0.  ,  0.  ,  0.  ],
                [ 0.  ,  0.  ,  0.  ,  0.  ,  1.7 ,  0.  ,  0.  ,  0.  ],
                [ 0.  ,  0.  ,  0.  ,  0.  ,  1.55,  0.  ,  0.  ,  0.  ],
                [ 0.  ,  1.7 ,  0.  ,  0.  ,  1.7 ,  0.  ,  0.  ,  0.  ],
                [ 0.  ,  1.55,  0.  ,  1.55,  1.55,  0.  ,  0.  ,  0.  ],
                [ 0.  ,  1.7 ,  0.  ,  0.  ,  1.7 ,  0.  ,  0.  ,  0.  ],
                [ 0.  ,  0.  ,  0.  ,  1.55,  1.55,  0.  ,  0.  ,  0.  ],
                [ 0.  ,  1.55,  0.  ,  1.55,  1.55,  0.  ,  0.  ,  0.  ],
                [ 0.  ,  1.7 ,  0.  ,  0.  ,  1.7 ,  0.  ,  0.  ,  0.  ],
                [ 0.  ,  0.  ,  0.  ,  1.55,  1.55,  0.  ,  0.  ,  0.  ],
                [ 0.  ,  1.55,  0.  ,  1.55,  1.55,  0.  ,  0.  ,  0.  ],
                [ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  2.  ,  0.  ,  0.  ],
                [ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  2.  ,  0.  ],
                [ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  2.  ],
                [ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  2.  ],
                [ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  2.  ],
                [ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  2.  ],
                [ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  2.  ],
                [ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  2.  ]]))
        """

        mol_channel = channel_init[:, [CHANNELS_INIT["Hydrophobic"],
                                       CHANNELS_INIT["Aromatic"],
                                       CHANNELS_INIT["Acceptor"],
                                       CHANNELS_INIT["Donor"],
                                       CHANNELS_INIT["HeavyAtom"]]]
        n_atoms = len(mol_coords)
        # treat aromatic/acceptor/donor as atoms, channel radii = PHARM_RADIUS
        aromatic_centers = self._get_aromatic_centers(
            mol_coords,
            self._get_aromatic_groups(rdkit_mol))
        if IS_SHUFFLE:
            aromatic_centers = self._shuffle_loc(aromatic_centers)
        n_aromatic = len(aromatic_centers)

        acceptor_rows = (mol_channel[:, CHANNELS["Acceptor"]] > 0.01)
        acceptor_coords = mol_coords[acceptor_rows]
        if IS_SHUFFLE:
            acceptor_coords = self._shuffle_loc(acceptor_coords)
        n_acceptor = len(acceptor_coords)

        donor_rows = (mol_channel[:, CHANNELS["Donor"]] > 0.01)
        donor_coords = mol_coords[donor_rows]
        if IS_SHUFFLE:
            donor_coords = self._shuffle_loc(donor_coords)

        # append phamacophore "atoms" to the molecular coordinates and channels
        coords_all = np.vstack([mol_coords, aromatic_centers, acceptor_coords,
                               donor_coords])

        channel_all = np.zeros((len(coords_all), NCHANNELS))
        channel_all[:n_atoms, :mol_channel.shape[1]] = mol_channel
        channel_all[n_atoms:(n_atoms + n_aromatic),
                    CHANNELS["Pharm_Aromatic"]] = PHARM_RADIUS
        channel_all[(n_atoms + n_aromatic):(n_atoms + n_aromatic + n_acceptor),
                    CHANNELS["Pharm_Acceptor"]] = PHARM_RADIUS
        channel_all[(n_atoms + n_aromatic + n_acceptor):,
                    CHANNELS["Pharm_Donor"]] = PHARM_RADIUS
        return coords_all, channel_all

    def voxelize(self):
        """
        This function generates a voxelized representation of a conformer
        Each grid has one value for each channel

        Returns
        -------------
        voxels : 4D numpy array [NVOXEL_X][NVOXEL_Y][NVOXEL_Z][NCHANNELS]
            store the voxel value for each property

        Other Parameters
        -------------
        voxel_values : 2D numpy.array
            shape(NVOXEL_X* NVOXEL_Y* NVOXEL_Z, NCHANNELS)
        voxel_centers : 2D numpy.array, shape(NVOXEL_X* NVOXEL_Y* NVOXEL_Z, 3)
            the "center" of each grid
            - the coordinates of the bottom-left corner
        nvoxels : numpy.ndarray, [NVOXEL_X, NVOXEL_Y, NVOXEL_Z]

        Examples
        -------------
        >>> VoxelMol('c1ccccc1').voxelize().shape
        (32, 32, 32, 8)

        See also
        -------------
        revoxelize
        """

        rdkit_mol = self.mol
        smallmol = toSmallMol(rdkit_mol)
        center_mol = get_center(smallmol)
        coords_mol = get_coords(smallmol)

        # generate coordinates and channels with pharmacophores added
        channel_mol = mk_voxel.getChannels(smallmol)[0]
        if ADD_PHARM:
            coords_all, channel_all = self.add_pharmacophore(
                coords_mol, channel_mol, rdkit_mol)
        else:
            coords_all, channel_all = coords_mol, channel_mol
        # print(self.add_pharmacophore(coords_mol, channel_mol, rdkit_mol))

        # do random rotation and 2A translation, including pharmacophores
        if IS_SHUFFLE:
            rot_matrix = uniformRandomRotation()
            coords_all = np.dot((coords_all - center_mol),
                                np.transpose(rot_matrix)) + center_mol
            center_mol = self._shuffle_loc(center_mol, by_range=DISPLACEMENT)

        # generate voxel descriptors (NCHANNELS channels) of each grid
        # for randomly transformed coordinates inside a pre-defined box
        voxel_centers, nvoxels = mk_voxel.getCenters(
            mol=smallmol,
            boxsize=[NVOXEL_X * VOXELSIZE,
                     NVOXEL_Y * VOXELSIZE,
                     NVOXEL_Z * VOXELSIZE],
            center=center_mol,
            voxelsize=VOXELSIZE)
        ncenters = voxel_centers.shape[0]
        natoms = coords_all.shape[0]
        nchannels = channel_all.shape[1]
        trunc = TRUNC * TRUNC

        pymp.config.nested = True
        voxel_values = pymp.shared.array((ncenters, nchannels))
        with pymp.Parallel(N_PROC) as p1:
            # for atom_row in range(natoms):
            for atom_row in p1.range(0, natoms):
                atom = coords_all[atom_row, :]
                atomradii = channel_all[atom_row, :]
                for grid_row in range(ncenters):
                    grid = voxel_centers[grid_row, :]
                    dist = np.zeros(3)
                    dist[0] = atom[0] - grid[0]
                    dist[1] = atom[1] - grid[1]
                    dist[2] = atom[2] - grid[2]
                    dist2 = (dist[0] * dist[0] +
                             dist[1] * dist[1] +
                             dist[2] * dist[2])
                    if dist2 > trunc:
                        continue
                    dist_r = 1 / sqrt(dist2)
                    for channel_col in range(nchannels):
                        if atomradii[channel_col] == 0:
                            continue
                        x = atomradii[channel_col] * dist_r
                        x3 = x * x * x
                        x12 = x3 * x3 * x3 * x3
                        value = 1 - exp(-x12)
                        voxel_values[grid_row, channel_col] = max(
                            voxel_values[grid_row, channel_col],
                            value)
        return voxel_values.reshape(NVOXEL_X, NVOXEL_Y, NVOXEL_Z,
                                    -1).astype(np.float32)

    def revoxelize(self):
        self.voxels = self.voxelize()

    def view_voxels(self, voxelsize=VOXELSIZE, draw='solid'):
        """
        This function calls VMD to visualize the voxel features,
        adapted from mk_voxel.viewVoxelFeatures

        Parameters
        -------------
        voxelsize : int
        draw : str
            vmd visualization method: 'solid', 'wireframe', 'points'
        """
        from moleculekit.vmdgraphics import VMDIsosurface

        voxelsize = np.repeat(voxelsize, 3)
        features = self.voxels
        smallmol = toSmallMol(self.mol)
        center_mol = get_center(smallmol)
        voxel_centers, nvoxels = mk_voxel.getCenters(
            mol=smallmol,
            boxsize=[NVOXEL_X * VOXELSIZE,
                     NVOXEL_Y * VOXELSIZE,
                     NVOXEL_Z * VOXELSIZE],
            center=center_mol,
            voxelsize=VOXELSIZE)
        centers = voxel_centers.reshape(list(nvoxels) + [3, ])
        loweredge = np.min(centers, axis=(0, 1, 2)) - (voxelsize / 2)

        channel_dict = CHANNELS if ADD_PHARM else CHANNELS_INIT
        for i, channel in enumerate(channel_dict):
            VMDIsosurface(features[..., i], loweredge, voxelsize, isovalue=0.5,
                          color=i, name=channel, draw=draw)

    def get_index_by_channel(self, channel):
        """
        This function returns an array of indices of the grids which feature
        the spefified channel

        Parameters
        -------------
        channel : str or int
            if str, it must be one of the following:
                CHANNELS = {"Hydrophobic": 0, "Aromatic": 1, "Acceptor": 2,
                            "Donor": 3, "HeavyAtom": 4, "Pharm_Aromatic": 5,
                            "Pharm_Acceptor": 6, "Pharm_Donor": 7}
            if int, it should be within 0 to 7

        Returns
        -------------
        indices : 2D numpy array of shape(#grids, 3)

        Examples
        -------------
        >>> VoxelMol('c1ccccc1').get_index_by_channel(0).shape[0]
        898
        """
        channel_dict = CHANNELS if ADD_PHARM else CHANNELS_INIT
        if isinstance(channel, str):
            try:
                channel = channel_dict[channel]
            except KeyError:
                print('Channel does not exist.')
                return None
        if not isinstance(channel, int) or channel >= len(channel_dict):
            print('Channel does not exist.')
            return None
        indices = []
        voxels = self.voxels
        for x in range(voxels.shape[0]):
            for y in range(voxels.shape[1]):
                for z in range(voxels.shape[2]):
                    if voxels[x][y][z][channel] > 0:
                        indices.append([x, y, z])
        return np.asarray(indices)

    def num_grids_per_channel(self):
        """
        This function returns the number of grids that feature each channel

        Returns
        -------------
        grid_dict : dict

        Examples
        -------------
        >>> VoxelMol('c1ccccc1').num_grids_per_channel()
        {'Hydrophobic': 904,
         'Aromatic': 904,
         'Acceptor': 0,
         'Donor': 0,
         'HeavyAtom': 904,
         'Pharm_Aromatic': 514,
         'Pharm_Acceptor': 0,
         'Pharm_Donor': 0}
        """

        grid_dict = {}
        channel_dict = CHANNELS if ADD_PHARM else CHANNELS_INIT
        for i in range(NCHANNELS):
            grid_dict[list(channel_dict.keys())[i]] = (
                self.get_index_by_channel(i).shape[0])
        return grid_dict


def toSmallMol(mol):
    """
    This is a wrapper function that converts a rdkit mol object to
    moleculekit SmallMol without adding H by default

    Parameters
    -------------
    mol : rdkit.Chem.rdchem.Mol object with one conformation
        or VoxelMol object

    Returns
    -------------
    mk_mol : moleculekit.smallmol.smallmol.SmallMol object
    """

    if isinstance(mol, VoxelMol):
        mol = mol.mol
    if isinstance(mol, Chem.Mol):
        return mk.SmallMol(mol, fixHs=False)
    return None


def get_coords(mol):
    """
    This function returns the coordinates of each atom in the molecule

    Parameters
    -------------
    mol : rdkit.Chem.rdchem.Mol object
        or moleculekit.smallmol.smallmol.SmallMol object
        with one conformation

    Returns
    -------------
    coords : 2D numpy array of shape(#atoms, 3)

    Examples
    -------------
    >>> get_coords(VoxelMol('c1ccccc1').mol)
    array([[ 0.80675215, -1.13780083,  0.00997585],
           [ 1.38863303,  0.12981945,  0.01964016],
           [ 0.58188063,  1.26762022,  0.00966327],
           [-0.80675193,  1.13780056, -0.00997713],
           [-1.38863321, -0.12981974, -0.01964149],
           [-0.58188072, -1.26762027, -0.00966646]])
    """

    if isinstance(mol, mk.SmallMol):
        return mol.get('coords')[:, :, 0]
    elif isinstance(mol, Chem.Mol):
        return mol.GetConformer().GetPositions()
    else:
        return np.empty(shape=(0, 3))


def get_center(mol):
    """
    This function returns the coordinates of the center of the molecule

    Parameters
    -------------
    mol : rdkit.Chem.rdchem.Mol object
        or moleculekit.smallmol.smallmol.SmallMol object
        with one conformation

    Returns
    -------------
    center : 1D numpy array of [x, y, z]

    Examples
    -------------
    >>> get_center(VoxelMol('c1ccccc1').mol)
    array([ -1.09883585e-08,  -1.00937274e-07,  -9.66406850e-07])
    """

    if isinstance(mol, mk.SmallMol):
        return mol.getCenter()
    elif isinstance(mol, Chem.Mol):
        points = Chem.rdMolTransforms.ComputeCentroid(mol.GetConformer())
        return np.array([points.x, points.y, points.z])
    else:
        return np.array([])
