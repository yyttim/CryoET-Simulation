import numpy as np
import random

class SpecialDistributions:
    @staticmethod
    def gaussian_distribution(mean, std_dev, size=1):
        return np.random.normal(mean, std_dev, size)

    @staticmethod
    def uniform_distribution(lower, upper, size=1):
        return np.random.uniform(lower, upper, size)

    @staticmethod
    def exponential_distribution(scale, size=1):
        return np.random.exponential(scale, size)

    @staticmethod
    def custom_distribution(custom_func, size=1):
        return custom_func(size)

class OccGen():
    def __init__(self, occ_rg):
        assert hasattr(occ_rg, '__len__') and (len(occ_rg) == 2) and (occ_rg[0] <= occ_rg[1])
        assert (occ_rg[0] >= 0) and (occ_rg[1] <= 100)
        self.__occ_rg = occ_rg

    def gen_occupancy(self):
        return random.uniform(self.__occ_rg[0], self.__occ_rg[1])

    def insert_molecules(self, molecule_list, distribution_method, **kwargs):
        insertions = []
        for molecule in molecule_list:
            if distribution_method == 'gaussian':
                pos = SpecialDistributions.gaussian_distribution(kwargs['mean'], kwargs['std_dev'])
            elif distribution_method == 'uniform':
                pos = SpecialDistributions.uniform_distribution(kwargs['lower'], kwargs['upper'])
            elif distribution_method == 'exponential':
                pos = SpecialDistributions.exponential_distribution(kwargs['scale'])
            elif distribution_method == 'custom':
                pos = SpecialDistributions.custom_distribution(kwargs['custom_func'])
            else:
                pos = self.gen_occupancy()
            insertions.append((molecule, pos))
        return insertions
