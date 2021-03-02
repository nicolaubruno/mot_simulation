from distutils.core import setup, Extension

mot_sim_module = Extension(
    'mot_sim',\
    sources = ['model/C_extension/mot_sim.c', 'model/C_extension/Py_wrapper.c']
)

setup(
    name = 'mot_sim',\
    version = '2.0',\
    description = 'C Extension with optimized functions related to the MOTSim',\
    ext_modules = [mot_sim_module]
)