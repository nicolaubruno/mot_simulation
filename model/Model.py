#
# Libraries and modules
import pandas as pd

#
class Model:

    ''' Attributes '''

    # DataFrame with all atoms
    @property
    def all_atoms(self):
        return self._all_atoms

    # Index of the attribute all_atoms corresponding to the atom that will be used in the simulation (active atom)
    @property
    def active_atom_id(self):
        return self._active_atom_id

    # Planck's constant
    @property
    def h(self):
        return 6.62607004e-34 # J * s

    # Elementary charge
    @property
    def e(self):
        return 1.60217662e-19 # C

    # Speed of light
    @property
    def c(self):
        return 2.99792458e8 # m / s

    # Boltzmann constant
    @property
    def k_B(self):
        return 1.38064852e-23 # J / K

    # Bohr magneton
    @property
    def mu_B(self):
        return 9.274009994e-24 # J / T
    

    ''' Methods '''

    #
    def __init__(self):
        atoms_data = {
            'symbol': ['Dy', 'Er'],\
            'name': ['Dysprosium', 'Erbium'],\
            'atomic_number': [66, 68],\
            'mass': [163.92917475, 166.0]
        }

        self._all_atoms = pd.DataFrame(atoms_data, index=atoms_data['symbol'])

        print(self._all_atoms)