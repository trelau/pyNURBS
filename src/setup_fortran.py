import setuptools
from ConfigParser import ConfigParser
from numpy.distutils.core import setup, Extension

# Read configuration file and set defaults
config = ConfigParser()
config.add_section('extra_f90_compile_args')
config.set('extra_f90_compile_args', 'optimize', '-O0')
config.set('extra_f90_compile_args', 'debug', '0')

config.add_section('extra_f77_compile_args')
config.set('extra_f77_compile_args', 'optimize', '-O0')
config.set('extra_f77_compile_args', 'debug', '0')
config.read('setup.cfg')

# Config ,utilities, and optimization
util_src = ['src/config.f90',
            'src/utils/compare.f90',
            'src/utils/math.f90',
            'src/utils/generic_list.f90',
            'src/opt/nelder_mead.f90']
util_only = ['dummy_compare', 'dummy_math', 'dummy_list', 'dummy_config',
             'dummy_simplex']

# MINPACK
minpack_src = ['src/minpack/dogleg.f',
               'src/minpack/dpmpar.f',
               'src/minpack/enorm.f',
               'src/minpack/fdjac1.f',
               'src/minpack/qform.f',
               'src/minpack/qrfac.f',
               'src/minpack/r1mpyq.f',
               'src/minpack/r1updt.f',
               'src/minpack/hybrd.f',
               'src/minpack/hybrd1.f',
               'src/minpack/hybrj.f',
               'src/minpack/hybrj1.f']

# Geometry
geom_src = ['src/geometry/geom_utils.f90',
            'src/geometry/geom_results.f90',
            'src/geometry/evaluate.f90',
            'src/geometry/modify.f90',
            'src/geometry/divide.f90',
            'src/geometry/calculate.f90',
            'src/geometry/project.f90',
            'src/geometry/invert.f90',
            'src/geometry/bounding_box.f90',
            'src/geometry/intersect_bbox.f90',
            'src/geometry/intersect_curve.f90',
            'src/geometry/intersect_triangle.f90',
            'src/geometry/tessellate.f90',
            'src/geometry/intersect_surface.f90',
            'src/geometry/map_surface.f90']
geom_only = ['dehomogenize_array1d', 'dehomogenize_array2d',
             'equivalence_array1d',
             'curve_point', 'curve_points', 'bezier_curve_points',
             'surface_point', 'surface_points', 'bezier_surface_points',
             'rat_curve_derivs', 'rat_surface_derivs',
             'decasteljau1', 'decasteljau2', 'basis_funs',
             'ders_basis_funs', 'find_span', 'curve_derivs_alg1',
             'surface_derivs_alg1', 'curve_knot_ins', 'surface_knot_ins',
             'decompose_curve', 'decompose_surface', 'is_curve_flat',
             'is_surface_flat', 'split_bezier_curve',
             'extract_bezier_curve', 'split_bezier_surface',
             'cp_lengths', 'cp_net_areas', 'arc_length_bezier',
             'surface_area_bezier',
             'project_point_to_curve', 'project_point_to_surface',
             'dummy_bbox', 'dummy_bi', 'intersect_curve_curve',
             'intersect_curve_plane', 'intersect_curve_surface',
             'intersect_triangle_ray', 'intersect_triangles',
             'intersect_triangle_plane', 'adaptive_curve_tessellate',
             'adaptive_surface_tessellate', 'intersect_surface_plane',
             'intersect_surface_surface',
             'refine_spi_point', 'refine_ssi_point', 'icurve_eval',
             'find_mult_knots', 'invert_points_on_plane',
             'surface_distance_map']

# Fortran flags
f90_opt = config.get('extra_f90_compile_args', 'optimize')
f90_flags = ['-Wall', '-Wextra', '-Wconversion', '-pedantic',
             '-Wunderflow']

if config.get('extra_f90_compile_args', 'debug') == '1':
    f90_debug = ['-g', '-fbacktrace', '-fcheck=all',
                 '-ffpe-trap=invalid,zero,overflow,underflow',
                 '-ffpe-summary=all']
else:
    f90_debug = []
f90_args = [f90_opt] + f90_flags + f90_debug

f77_opt = config.get('extra_f77_compile_args', 'optimize')
f77_flags = ['-Wall', '-Wextra']
if config.get('extra_f77_compile_args', 'debug') == '1':
    f77_debug = ['-g', '-fbacktrace', '-fcheck=all',
                 '-ffpe-trap=invalid,zero,overflow,underflow',
                 '-ffpe-summary=all']
else:
    f77_debug = []
f77_args = [f77_opt] + f77_flags + f77_debug

# Build _lib extension
src = util_src + minpack_src + geom_src
only = util_only + geom_only
_ext = Extension('pynurbs.lib.pynurbs', src,
                 f2py_options=['only:'] + only + [':'],
                 extra_f90_compile_args=f90_args,
                 extra_f77_compile_args=f77_args)

kwds = {'name': 'pyNURBS',
        'version': '1.0',
        'package_dir': {'pyNURBS': 'pynurbs'},
        'packages': setuptools.find_packages('.'),
        'install_requires': ['numpy', 'scipy', 'mayavi'],
        'author': 'trelau',
        'description': 'Simple Python-based NURBS library.',
        'license': 'BSD',
        'ext_modules': [_ext]}

setup(**kwds)
