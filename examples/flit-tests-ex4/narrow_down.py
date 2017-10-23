#!/usr/bin/env python3

import sys
from glob import glob
import multiprocessing as mp
import subprocess
import tempfile
import itertools as it

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

src_of_interest = [
    '../../fem/fe.cpp',
    '../../fem/intrules.cpp',
    '../../linalg/densemat.cpp',
    '../../mesh/mesh.cpp',
    '../../mesh/nurbs.cpp',
    'tests/Example04.cpp',
    ]

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

def only_one_dev():
    'One file in executable, all others in shared library'
    dev_src_list = []
    gt_src_list = []
    for src_file in src:
        dev_src_list.append(set([src_file]))
        gt_src_list.append(set(src) - dev_src_list[-1])
    print('Running tests')
    with mp.Pool() as p:
        results = p.starmap(run_test, zip(dev_src_list, gt_src_list))
    for dev_src, success in zip(dev_src_list, results):
        if success:
            print('\nCompiling all but {0} into libgt.so'.format(dev_src) +
                  '\n  Caused a success')
        else:
            print('.', end='')
    print('\nFully finished')

def all_two_combos_dev():
    'Two files in executable, all others in shared library'
    dev_src_list = []
    gt_src_list = []
    already_paired = set()
    for src_file in src:
        for src_file_2 in src:
            if src_file == src_file_2 \
               or (src_file, src_file_2) in already_paired \
               or (src_file_2, src_file) in already_paired:
                continue
            dev_src_list.append(set([src_file, src_file_2]))
            gt_src_list.append(set(src) - dev_src_list[-1])
            already_paired.add((src_file, src_file_2))
    print('Running tests')
    with mp.Pool() as p:
        results = p.starmap(run_test, zip(dev_src_list, gt_src_list))
    for dev_src, success in zip(dev_src_list, results):
        if success:
            print('\nCompiling all but {0} into libgt.so'.format(dev_src) +
                  '\n  Caused a success')
        else:
            print('.', end='')
    print('\nFully finished')

def only_one_gt():
    'One file in shared library, all others in executable'
    dev_src_list = []
    gt_src_list = []
    for src_file in src:
        dev_src_list.append(set(src) - set([src_file]))
        gt_src_list.append(set(src) - dev_src_list[-1])
    print('Running tests')
    with mp.Pool() as p:
        results = p.starmap(run_test, zip(dev_src_list, gt_src_list))
    for gt_src, success in zip(gt_src_list, results):
        if not success:
            print('\nCompiling {0} into libgt.so'.format(gt_src) +
                  '\n  Caused a failure')
        else:
            print('.', end='')
    print('\nFully finished')

def all_interesting_files():
    '''
    Combos of putting interesting files in shared library, all others in
    executable
    '''
    print('Interesting files investigated:')
    for filename in src_of_interest:
        print(' ', filename)
    def powerset(A):
        iterable = it.chain.from_iterable(
                it.combinations(A, r) for r in range(len(A) + 1))
        for elem in iterable:
            yield set(elem)
    combos = powerset(src_of_interest)
    dev_src_list = []
    gt_src_list = []
    src_set = set(src)
    for gt in combos:
        gt_src_list.append(gt)
        dev_src_list.append(src_set - gt)
    with mp.Pool() as p:
        results = p.starmap(run_test, zip(dev_src_list, gt_src_list))
    for gt_src, success in zip(gt_src_list, results):
        if not success:
            print('\nCompiling {0} into libgt.so'.format(gt_src) +
                  '\n  Caused a failure')
        else:
            print('.', end='')
    print('\nFully finished')

def main(arguments):
    'Main entry point'
    print('Compiling objects before experiments')
    compile_objects(src)
    print('Finished compiling objects')
    print()

    #print('Running only_one_gt()')
    #only_one_gt()
    #print()

    #print('Running only_one_dev()')
    #only_one_dev()
    #print()

    #print('Running all_two_combos_dev()')
    #all_two_combos_dev()
    #print()

    print('Running all_interesting_files()')
    all_interesting_files()
    print()

if __name__ == '__main__':
    main(sys.argv[1:])
