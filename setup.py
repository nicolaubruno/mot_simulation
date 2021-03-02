from distutils.core import setup, Extension

sim_atom_module = Extension(
    'sim_atom',\
    sources = ['model/C_extensions/wrapper.c', 'model/C_extensions/sim_atom.c']
)

setup(
    name = 'sim_atom',\
    version = '2.0',\
    description = 'Simulation of a single atom',\
    ext_modules = [sim_atom_module]
)