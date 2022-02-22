"""
Classes for modeling networks polymers.
A networks is a combination of a polymer in a volume
"""

__author__ = 'Antonio Martinez-Sanchez'

import scipy as sp

from polnet.poly import *
from polnet.polymer import SAWLC, SAWLCPoly, HelixFiber
from polnet.lrandom import PGen, PGenHelixFiber, PGenHelixFiberB
from abc import ABC, abstractmethod


class Network(ABC):
    """
    General class for a network of polymers
    """

    def __init__(self, voi, v_size, svol=None):
        """
        Construction
        :param voi: a 3D numpy array to define a VOI (Volume Of Interest) for polymers
        :param v_size: voxel size (default 1)
        :param svol: monomer subvolume as a numpy ndarray (default None)
        """
        self.set_voi(voi)
        self.__vol = (self.__voi > 0).sum() * v_size * v_size * v_size
        self.__v_size = v_size
        self.__pl_occ = 0
        self.__pl = list()
        self.__svol = svol
        if self.__svol is not None:
            assert isinstance(self.__svol, np.ndarray)

    def get_pmers_list(self):
        return self.__pl

    def get_polymer_occupancy(self):
        return self.__pl_occ

    def add_polymer(self, polymer):
        self.__pl.append(polymer)
        self.__pl_occ += 100. * (polymer.get_vol() / self.__vol)
        print('Occ: ', self.__pl_occ)

    @abstractmethod
    def build_network(self):
        """
        Builds an instance of the network
        :return: None        """
        raise NotImplemented

    def get_voi(self):
        """
        Get the VOI
        :return: an ndarray
        """
        return self.__voi

    def get_gtruth(self, thick=1):
        """
        Get the ground truth tomogram
        :param thick: ground truth tickness in voxels (default 1)
        :return: a binary numpy 3D array
        """
        hold_gtruth = self.gen_vtp_points_tomo()
        if thick >= 1:
            hold_gtruth = sp.ndimage.morphology.binary_dilation(hold_gtruth, iterations=int(thick))
        return hold_gtruth

    def set_voi(self, voi):
        """
        Set the VOI
        :param voi:
        """
        assert isinstance(voi, np.ndarray)
        if voi.dtype is bool:
            self.__voi = voi
        else:
            self.__voi = voi > 0

    def get_vtp(self):
        """
        Get Polymers Network as a vtkPolyData with their surfaces
        :return: a vtkPolyData
        """

        app_flt = vtk.vtkAppendPolyData()

        # Polymers loop
        for pol in self.__pl:
            app_flt.AddInputData(pol.get_vtp())
        app_flt.Update()

        return app_flt.GetOutput()

    def get_skel(self, add_verts=True, add_lines=True, verts_rad=0):
        """
        Get the polymer as a skeleton, each monomer is a point or sphere and lines connecting monomers
        :param add_verts: if True (default) the vertices are included in the vtkPolyData
        :param add_lines: if True (default) the lines are included in the vtkPolyData
        :param verts_rad: if verts is True then sets the vertex radius, if <=0 a vertices are just points
        :return: a vtkPolyData
        """
        app_flt = vtk.vtkAppendPolyData()

        # Polymers loop
        for pol in self.__pl:
            app_flt.AddInputData(pol.get_skel(add_verts, add_lines, verts_rad))

        # Update and return
        app_flt.Update()
        return app_flt.GetOutput()

    def gen_vtp_points_tomo(self):
        """
        Generates a binary tomogram where True elements correspond with the polydata closes voxel projection
        :return: a binary VOI shaped numpy array
        """
        nx, ny, nz = self.__voi.shape
        hold_tomo = np.zeros(shape=(nx, ny, nz), dtype=bool)
        hold_vtp_skel = self.get_skel()
        for i in range(hold_vtp_skel.GetNumberOfPoints()):
            x, y, z = hold_vtp_skel.GetPoint(i)
            x, y, z = int(round(x)), int(round(y)), int(round(z))
            if (x >= 0) and (y >= 0) and (z >= 0) and  (x < nx) and (y < ny) and (z < nz):
                hold_tomo[x, y, z] = True
        return hold_tomo

    def insert_density_svol(self, m_svol, tomo, v_size=1, merge='max', off_svol=None):
        """
        Insert a polymer network as set of subvolumes into a tomogram
        :param m_svol: input monomer sub-volume reference
        :param tomo: tomogram where m_svol is added
        :param v_size: tomogram voxel size (default 1)
        :param merge: merging mode, valid: 'min' (default), 'max', 'sum' and 'insert'
        :param off_svol: offset coordinates for sub-volume monomer center coordinates
        """
        assert isinstance(m_svol, np.ndarray) and (len(m_svol.shape) == 3)
        assert isinstance(tomo, np.ndarray) and (len(tomo.shape) == 3)
        assert (merge == 'max') or (merge == 'min') or (merge == 'sum') or (merge == 'insert')
        assert v_size > 0
        if off_svol is not None:
            assert isinstance(off_svol, np.ndarray) and (len(off_svol) == 3)

        for pl in self.__pl:
            pl.insert_density_svol(m_svol, tomo, v_size, merge=merge, off_svol=off_svol)

    def add_monomer_to_voi(self, mmer, mmer_svol=None):
        """
        Adds a monomer to VOI mask
        :param mmer: monomer to define rigid transformations
        :param mmer_voi: subvolume (binary numpy ndarray) with monomer VOI
        """
        assert isinstance(mmer_svol, np.ndarray) and (mmer_svol.dtype == bool)
        v_size_i = 1. / self.__v_size
        tot_v = np.asarray((0., 0., 0.))
        hold_svol = mmer_svol > 0
        for trans in mmer.get_trans_list():
            if trans[0] == 't':
                tot_v += (trans[1] * v_size_i)
            elif trans[0] == 'r':
                hold_svol = tomo_rotate(hold_svol, trans[1], order=0, mode='constant', cval=hold_svol.max(),
                                        prefilter=False)
                # hold_svol = tomo_rotate(hold_svol, trans[1], mode='constant', cval=hold_svol.min())
        insert_svol_tomo(hold_svol, self.__voi, tot_v, merge='min')


class NetSAWLC(Network):
    """
    Class for generating a network of SAWLC polymers
    """

    def __init__(self, voi, v_size, l_length, m_surf, max_p_length, gen_pol_lengths, occ, over_tolerance=0,
                 poly=None, svol=None):
        """
        Construction
        :param voi: a 3D numpy array to define a VOI (Volume Of Interest) for polymers
        :param v_size: voxel size (default 1)
        :param l_length: polymer link length
        :param m_surf: monomer surf
        :param max_p_length: maximum polymer length
        :param gen_pol_lengths: a instance of a random generation model (random.PGen) to determine the achievable
        lengths for polymers
        :param occ: occupancy threshold in percentage [0, 100]%
        :param over_tolerance: fraction of overlapping tolerance for self avoiding (default 0, in range [0,1))
        :param poly: it allows to restrict monomer localizations to a polydata
        :param svol: monomer subvolume as a numpy ndarray (default None)
        :param off_svol: offset coordinates in voxels for shifting sub-volume monomer center coordinates (default None)
        """

        # Initialize abstract varibles
        super(NetSAWLC, self).__init__(voi, v_size, svol=svol)

        # Input parsing
        assert l_length > 0
        assert isinstance(m_surf, vtk.vtkPolyData)
        assert max_p_length >= 0
        assert isinstance(gen_pol_lengths, PGenHelixFiber)
        assert (occ >= 0) and (occ <= 100)
        assert (over_tolerance >= 0) and (over_tolerance <= 100)
        if poly is not None:
            assert isinstance(poly, vtk.vtkPolyData)

        # Variables assignment
        self.__l_length, self.__m_surf = l_length, m_surf
        self.__max_p_length = max_p_length
        self.__gen_pol_lengths = gen_pol_lengths
        self.__occ, self.__over_tolerance = occ, over_tolerance
        self.__poly = poly

    def build_network(self):
        """
        Add polymers following SAWLC model until an occupancy limit is passed
        :return:
        """

        # Computes the maximum number of tries
        c_try = 0
        n_tries = math.ceil(self._Network__vol / poly_volume(self.__m_surf))

        # Network loop
        while (c_try <= n_tries) and (self._Network__pl_occ < self.__occ):

            # Polymer initialization
            c_try += 1
            p0 = np.asarray((self._Network__voi.shape[0] * self._Network__v_size * random.random(),
                             self._Network__voi.shape[1] * self._Network__v_size * random.random(),
                             self._Network__voi.shape[2] * self._Network__v_size * random.random()))
            max_length = self.__gen_pol_lengths.gen_length(0, self.__max_p_length)
            if self.__poly is None:
                hold_polymer = SAWLC(self.__l_length, self.__m_surf, p0)
            else:
                hold_polymer = SAWLCPoly(self.__poly, self.__l_length, self.__m_surf, p0)
            if hold_polymer.get_monomer(-1).overlap_voi(self._Network__voi, self._Network__v_size):
                continue
            self.add_monomer_to_voi(hold_polymer.get_monomer(-1), self._Network__svol)

            # Polymer loop
            not_finished = True
            while (hold_polymer.get_total_len() < max_length) and not_finished:
                monomer_data = hold_polymer.gen_new_monomer(self.__over_tolerance, self._Network__voi,
                                                            self._Network__v_size)
                if monomer_data is None:
                    not_finished = False
                else:
                    new_len = points_distance(monomer_data[0], hold_polymer.get_tail_point())
                    if hold_polymer.get_total_len() + new_len < max_length:
                        # ) and (monomer_data[3].overlap_voi(self._Network__voi, self._Network__v_size)):
                        hold_polymer.add_monomer(monomer_data[0], monomer_data[1], monomer_data[2], monomer_data[3])
                        self.add_monomer_to_voi(hold_polymer.get_monomer(-1), self._Network__svol)
                    else:
                        not_finished = False

            # Updating polymer
            self.add_polymer(hold_polymer)
            # print('build_network: new polymer added with ' + str(hold_polymer.get_num_monomers()) +
            #       ' and length ' + str(hold_polymer.get_total_len()))


class NetHelixFiber(Network):
    """
    Class for generating a network of isolated helix fibers, unconnected and randomly distributed
    """

    def __init__(self, voi, v_size, l_length, m_surf, gen_hfib_params, occ, min_p_len, hp_len, mz_len, mz_len_f,
                 over_tolerance=0):
        """
        Construction
        :param voi: a 3D numpy array to define a VOI (Volume Of Interest) for polymers
        :param v_size: voxel size (default 1)
        :param l_length: polymer link length
        :param m_surf: monomer surf
        :param gen_hfib_params: a instance of a random generation model (random.PGen) to obtain random fiber
        parametrization
        :param min_p_len: minimum persistence length
        :param hp_len: helix period length
        :param mz_len: monomer length is z-axis
        :param mz_len_f: maximum length factor in z-axis
        :param occ: occupancy threshold in percentage [0, 100]%
        :param over_tolerance: fraction of overlapping tolerance for self avoiding (default 0, in range [0,1))
        """

        # Initialize abstract variables
        super(NetHelixFiber, self).__init__(voi, v_size)

        # Input parsing
        assert l_length > 0
        assert isinstance(m_surf, vtk.vtkPolyData)
        assert isinstance(gen_hfib_params, PGenHelixFiber)
        assert (occ >= 0) and (occ <= 100)
        assert (over_tolerance >= 0) and (over_tolerance <= 100)
        assert (min_p_len > 0) and (hp_len > 0) and (mz_len > 0) and (mz_len_f >= 0)

        # Variables assignment
        self.__l_length, self.__m_surf = l_length, m_surf
        self.__gen_hfib_params = gen_hfib_params
        self.__occ, self.__over_tolerance = occ, over_tolerance
        self.__min_p_len, self.__hp_len = min_p_len, hp_len
        self.__mz_len, self.__mz_len_f = mz_len, mz_len_f

    def build_network(self):
        """
        Add helix fibres until an occupancy limit is passed
        :return:
        """

        # Network loop
        while self._Network__pl_occ < self.__occ:

            # Polymer initialization
            p0 = np.asarray((self._Network__voi.shape[0] * self._Network__v_size * random.random(),
                             self._Network__voi.shape[1] * self._Network__v_size * random.random(),
                             self._Network__voi.shape[2] * self._Network__v_size * random.random()))
            max_length = math.sqrt(self._Network__voi.shape[0]**2 + self._Network__voi.shape[1]**2
                                  + self._Network__voi.shape[2]**2) * self._Network__v_size
            p_len = self.__gen_hfib_params.gen_persistence_length(self.__min_p_len)
            z_len_f = self.__gen_hfib_params.gen_zf_length(self.__mz_len_f)
            hold_polymer = HelixFiber(self.__l_length, self.__m_surf, p_len, self.__hp_len, self.__mz_len, z_len_f, p0)

            # Polymer loop
            not_finished = True
            while (hold_polymer.get_total_len() < max_length) and not_finished:
                monomer_data = hold_polymer.gen_new_monomer(self.__over_tolerance, self._Network__voi,
                                                            self._Network__v_size)
                if monomer_data is None:
                    not_finished = False
                else:
                    new_len = points_distance(monomer_data[0], hold_polymer.get_tail_point())
                    if hold_polymer.get_total_len() + new_len < max_length:
                        hold_polymer.add_monomer(monomer_data[0], monomer_data[1], monomer_data[2], monomer_data[3])
                    else:
                        not_finished = False

            # Updating polymer
            self.add_polymer(hold_polymer)
            # print('build_network: new polymer added with ' + str(hold_polymer.get_num_monomers()) +
            #       ' and length ' + str(hold_polymer.get_total_len()) + ': occupancy ' + str(self._Network__pl_occ))


class NetHelixFiberB(Network):
    """
    Class for generating a network of brancked helix fibers randomly distributed
    """

    def __init__(self, voi, v_size, l_length, m_surf, gen_hfib_params, occ, min_p_len, hp_len, mz_len, mz_len_f,
                 b_prop, max_p_branch=0, over_tolerance=0):
        """
        Construction
        :param voi: a 3D numpy array to define a VOI (Volume Of Interest) for polymers
        :param v_size: voxel size (default 1)
        :param l_length: polymer link length
        :param m_surf: monomer surf
        :param gen_hfib_params: a instance of a random generation model (random.PGen.NetHelixFiberB) to obtain random fiber
        parametrization for helix with branches
        :param occ: occupancy threshold in percentage [0, 100]%
        :param min_p_len: minimum persistence length
        :param hp_len: helix period length
        :param mz_len: monomer length is z-axis
        :param mz_len_f: maximum length factor in z-axis
        :param b_prob: branching probability, checked every time a new monomer is added
        :param max_p_branch: maximum number of branches per polymer, if 0 (default) then no branches are generated
        :param over_tolerance: fraction of overlapping tolerance for self avoiding (default 0, in range [0,1))
        """

        # Initialize abstract variables
        super(NetHelixFiberB, self).__init__(voi, v_size)

        # Input parsing
        assert l_length > 0
        assert isinstance(m_surf, vtk.vtkPolyData)
        assert isinstance(gen_hfib_params, PGenHelixFiberB)
        assert (occ >= 0) and (occ <= 100)
        assert (over_tolerance >= 0) and (over_tolerance <= 100)
        assert (min_p_len > 0) and (hp_len > 0) and (mz_len > 0) and (mz_len_f >= 0)
        assert (max_p_branch >= 0) and (b_prop >= 0)

        # Variables assignment
        self.__l_length, self.__m_surf = l_length, m_surf
        self.__gen_hfib_params = gen_hfib_params
        self.__occ, self.__over_tolerance = occ, over_tolerance
        self.__min_p_len, self.__hp_len = min_p_len, hp_len
        self.__mz_len, self.__mz_len_f = mz_len, mz_len_f
        self.__max_p_branch, self.__p_branches, self.__b_prop = max_p_branch, list(), b_prop

    def build_network(self):
        """
        Add helix fibres until an occupancy limit is passed
        :return:
        """

        # Network loop
        while self._Network__pl_occ < self.__occ:

            # Polymer initialization
            max_length = math.sqrt(self._Network__voi.shape[0]**2 + self._Network__voi.shape[1]**2
                                  + self._Network__voi.shape[2]**2) * self._Network__v_size
            p_len = self.__gen_hfib_params.gen_persistence_length(self.__min_p_len)
            z_len_f = self.__gen_hfib_params.gen_zf_length(self.__mz_len_f)
            branch = None
            if (self.__max_p_branch > 0) and self.__gen_hfib_params.gen_branch(self.__b_prop):
                branch = self.__gen_random_branch()
            if branch is None:
                p0 = np.asarray((self._Network__voi.shape[0] * self._Network__v_size * random.random(),
                                 self._Network__voi.shape[1] * self._Network__v_size * random.random(),
                                 self._Network__voi.shape[2] * self._Network__v_size * random.random()))
            else:
                p0 = branch.get_point()
            hold_polymer = HelixFiber(self.__l_length, self.__m_surf, p_len, self.__hp_len,
                                      self.__mz_len, z_len_f, p0)

            # Polymer loop
            not_finished = True
            while (hold_polymer.get_total_len() < max_length) and not_finished:
                monomer_data = hold_polymer.gen_new_monomer(self.__over_tolerance, self._Network__voi,
                                                            self._Network__v_size)
                if monomer_data is None:
                    not_finished = False
                else:
                    new_len = points_distance(monomer_data[0], hold_polymer.get_tail_point())
                    if hold_polymer.get_total_len() + new_len < max_length:
                        hold_polymer.add_monomer(monomer_data[0], monomer_data[1], monomer_data[2], monomer_data[3])
                    else:
                        not_finished = False

            # Updating polymer
            if branch is not None:
                self.add_polymer(hold_polymer)
                self.__p_branches.append(list())
                self.__add_branch(hold_polymer, branch)
            else:
                self.add_polymer(hold_polymer)
                self.__p_branches.append(list())
            # print('build_network: new polymer added with ' + str(hold_polymer.get_num_monomers()) +
            #       ' and length ' + str(hold_polymer.get_total_len()) + ': occupancy ' + str(self._Network__pl_occ))

    def get_branch_list(self):
        """
        Get all branches in a list
        :return: a single list with the branches
        """
        hold_list = list()
        for bl in self.__p_branches:
            for b in bl:
                hold_list.append(b)
        return hold_list

    def get_skel(self):
        """
        Get Polymers Network as a vtkPolyData as points and lines with branches
        :return: a vtkPolyData
        """

        # Initialization
        app_flt_l, app_flt_v, app_flt = vtk.vtkAppendPolyData(), vtk.vtkAppendPolyData(), vtk.vtkAppendPolyData()

        # Polymers loop
        p_type_l = vtk.vtkIntArray()
        p_type_l.SetName('type')
        p_type_l.SetNumberOfComponents(1)
        for pol in self._Network__pl:
            app_flt_l.AddInputData(pol.get_skel())
        app_flt_l.Update()
        out_vtp_l = app_flt_l.GetOutput()
        for i in range(out_vtp_l.GetNumberOfCells()):
            p_type_l.InsertNextTuple((1,))
        out_vtp_l.GetCellData().AddArray(p_type_l)

        # Branches loop
        p_type_v = vtk.vtkIntArray()
        p_type_v.SetName('type')
        p_type_v.SetNumberOfComponents(1)
        for i, branch in enumerate(self.get_branch_list()):
            app_flt_v.AddInputData(branch.get_vtp())
            # print('Point ' + str(i) + ': ' + str(branch.get_point()))
        app_flt_v.Update()
        out_vtp_v = app_flt_v.GetOutput()
        for i in range(out_vtp_v.GetNumberOfCells()):
            p_type_v.InsertNextTuple((2,))
        out_vtp_v.GetCellData().AddArray(p_type_v)

        # Merging branches and polymers
        app_flt.AddInputData(out_vtp_l)
        app_flt.AddInputData(out_vtp_v)
        app_flt.Update()

        return app_flt.GetOutput()

    def __gen_random_branch(self):
        """
        Generates a position point randomly for a branch on the filament network, no more than one branch per polymer
        :return: a branch
        """

        # Loop for polymers
        count, branch = 0, None
        while (count < len(self._Network__pl)) and (branch is None):
            hold_pid = random.randint(0, len(self._Network__pl)-1)
            if len(self.__p_branches[hold_pid]) < self.__max_p_branch:
                hold_pol = self._Network__pl[hold_pid]
                hold_mid = random.randint(0, len(hold_pol._Polymer__m) - 1)
                hold_m = hold_pol._Polymer__m[hold_mid]
                found = True
                for branch in self.__p_branches[hold_pid]:
                    if points_distance(hold_m.get_center_mass(), branch.get_point()) <= 2 * hold_m.get_diameter():
                        found = False
                if found:
                    branch = Branch(hold_m.get_center_mass(), hold_pid, hold_mid)
            count += 1

        return branch

    def __add_branch(self, polymer, branch):
        """
        Add a new branch to the polymer network
        :param polymer: targeting polymer where the branch is going to be added (starting polumer is obtained from
                        the branch)
        :param branch: branch to be added
        """
        branch.set_t_pmer(len(self._Network__pl) - 1)
        self.__p_branches[branch.get_s_pmer()].append(branch)


class Branch:
    """
    Class to model a brach in a Network
    """

    def __init__(self, point, s_pmer_id, s_mmer_id, t_pmer_id=None):
        """
        Constructor
        :param point: branch point coordinates (3-array)
        :param s_pmer_id: starting polymer ID
        :param s_mmer_id: starting monomer ID
        :param t_pmer_id: targeting polymer ID, (default None) it may be unknown at construction time, the targeting
                          monomer ID is 0.
        """
        assert hasattr(point, '__len__') and (len(point) == 3)
        assert (s_pmer_id >= 0) and (s_mmer_id >= 0)
        if t_pmer_id is not None:
            assert t_pmer_id >= 0
        self.__point = np.asarray(point, dtype=float)
        self.__s_pmer_id, self.__s_mmer_id, self.__t_pmer_id = s_pmer_id, s_mmer_id, t_pmer_id

    def set_t_pmer(self, t_pmer_id):
        """
        Set targeting polymer ID
        """
        assert t_pmer_id >= 0
        self.__t_pmer_id = t_pmer_id

    def get_point(self):
        """
        Get point coordinates
        """
        return self.__point

    def get_s_pmer(self):
        """
        Get starting polymer ID
        """
        return self.__s_pmer_id

    def get_s_mmer(self):
        """
        Get starting monomer ID
        """
        return self.__s_mmer_id

    def get_t_pmer(self):
        """
        Get targeting polymer ID
        """
        return self.__t_pmer_id

    def get_vtp(self):
        return point_to_poly(self.get_point())


