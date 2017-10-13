#!/usr/bin/env python3

import sys
from glob import glob
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

def run_test(dev_src, gt_src):
    '''
    Runs the two sets and returns True if the result is the same as the full
    dev build.
    '''
    infname = 'Makefile.in'
    with open(infname, 'r') as infile:
        content = str(infile.read())

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
            subprocess.check_call(['make', '-j30', '-f', outfname],
                                  stdout=subprocess.DEVNULL,
                                  stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError:
            print('Compiling {0} into libgt.so'.format(src_file) +
                  '\n  Caused a failure')
            return False
        else:
            return True

def main(arguments):
    'Main entry point'
    for src_file in src:
        dev_src = set(src) - set([src_file])
        gt_src = set(src) - dev_src
        print('Running test for {0}'.format(src_file))
        success = run_test(dev_src, gt_src)
        if not success:
            pass

if __name__ == '__main__':
    main(sys.argv[1:])
