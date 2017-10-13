#!/usr/bin/env python3

import sys
from glob import glob
import multiprocessing
import subprocess
import tempfile

src = [
    'main.cpp',
    '../../fem/bilinearform.cpp',
    '../../fem/bilininteg.cpp',
    '../../fem/coefficient.cpp',
    '../../fem/datacollection.cpp',
    '../../fem/eltrans.cpp',
    '../../fem/fe.cpp',
    '../../fem/fe_coll.cpp',
    '../../fem/fespace.cpp',
    '../../fem/geom.cpp',
    '../../fem/gridfunc.cpp',
    '../../fem/hybridization.cpp',
    '../../fem/intrules.cpp',
    '../../fem/linearform.cpp',
    '../../fem/lininteg.cpp',
    '../../fem/nonlinearform.cpp',
    '../../fem/nonlininteg.cpp',
    '../../fem/staticcond.cpp',
    '../../general/array.cpp',
    '../../general/error.cpp',
    '../../general/gzstream.cpp',
    '../../general/socketstream.cpp',
    '../../general/stable3d.cpp',
    '../../general/table.cpp',
    '../../general/tic_toc.cpp',
    '../../linalg/blockmatrix.cpp',
    '../../linalg/blockoperator.cpp',
    '../../linalg/blockvector.cpp',
    '../../linalg/densemat.cpp',
    '../../linalg/matrix.cpp',
    '../../linalg/ode.cpp',
    '../../linalg/operator.cpp',
    '../../linalg/solvers.cpp',
    '../../linalg/sparsemat.cpp',
    '../../linalg/sparsesmoothers.cpp',
    '../../linalg/vector.cpp',
    '../../mesh/element.cpp',
    '../../mesh/hexahedron.cpp',
    '../../mesh/mesh.cpp',
    '../../mesh/mesh_readers.cpp',
    '../../mesh/ncmesh.cpp',
    '../../mesh/nurbs.cpp',
    '../../mesh/point.cpp',
    '../../mesh/quadrilateral.cpp',
    '../../mesh/segment.cpp',
    '../../mesh/tetrahedron.cpp',
    '../../mesh/triangle.cpp',
    '../../mesh/vertex.cpp',
    ]

src.extend(glob('tests/*.cpp'))

def compile_objects(src):
    'Compiles shared targets such as object files and the devrun executable'
    assert(run_test(src, src, make_args=['-j30', 'objects'], quietout=False,
                    clean=False))
    assert(run_test(src, set(), make_args=['-j30', 'dev']))

def run_test(dev_src, gt_src, make_args=[], quietout=True, quieterr=True,
             clean=True):
    '''
    Runs the two sets and returns True if the result is the same as the full
    dev build.
    '''
    infname = 'Makefile.in'
    with open(infname, 'r') as infile:
        content = str(infile.read())

    extra_call_args = {}
    if quietout:
        extra_call_args['stdout'] = subprocess.DEVNULL
    if quieterr:
        extra_call_args['stderr'] = subprocess.DEVNULL

    replace_map = {}
    replace_map['DEV_SRC'] = '\n'.join(['DEV_SRC    += {0}'.format(x)
                                        for x in dev_src])
    replace_map['GT_SRC'] = '\n'.join(['GT_SRC     += {0}'.format(x)
                                       for x in gt_src])
    with tempfile.NamedTemporaryFile(dir=tempfile.gettempdir()) as outfile:
        outfname = outfile.name
        replace_map['Makefile'] = outfname
        outfile.write(bytes(content.format(**replace_map), 'utf-8'))
        try:
            subprocess.check_call(['make', '-f', outfname] + make_args,
                                  **extra_call_args)
        except subprocess.CalledProcessError:
            return False
        else:
            return True
        finally:
            if clean:
                subprocess.call(['make', 'clean', '-f', outfname],
                                **extra_call_args)

def main(arguments):
    'Main entry point'
    dev_src_list = []
    gt_src_list = []
    for src_file in src:
        dev_src_list.append(set(src) - set([src_file]))
        gt_src_list.append(set(src) - dev_src_list[-1])
    compile_objects(src)
    print('Finished compiling objects')
    print('Running tests')
    with multiprocessing.Pool() as p:
        results = p.starmap(run_test, zip(dev_src_list, gt_src_list))
    for gt_src, success in zip(gt_src_list, results):
        if not success:
            print('\nCompiling {0} into libgt.so'.format(gt_src) +
                  '\n  Caused a failure')
        else:
            print('.', end='')
    print('\nFully finished')

if __name__ == '__main__':
    main(sys.argv[1:])
