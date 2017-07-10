# This file is included at the end of the copied Makefile.  If you have some
# things you want to change about the Makefile, it is best to do it here.

# additional source files to compile other than what is in '.' and 'tests/'
# since those directories are added by a wildcard.
#SOURCE         += ../../fem/bilinearform.cpp
#SOURCE         += ../../fem/bilininteg.cpp
#SOURCE         += ../../fem/coefficient.cpp
#SOURCE         += ../../fem/datacollection.cpp
#SOURCE         += ../../fem/eltrans.cpp
#SOURCE         += ../../fem/fe.cpp
#SOURCE         += ../../fem/fe_coll.cpp
#SOURCE         += ../../fem/fespace.cpp
#SOURCE         += ../../fem/geom.cpp
#SOURCE         += ../../fem/gridfunc.cpp
#SOURCE         += ../../fem/hybridization.cpp
#SOURCE         += ../../fem/intrules.cpp
#SOURCE         += ../../fem/linearform.cpp
#SOURCE         += ../../fem/lininteg.cpp
#SOURCE         += ../../fem/nonlinearform.cpp
#SOURCE         += ../../fem/nonlininteg.cpp
#SOURCE         += ../../fem/staticcond.cpp        # 17/24 of the fem files
#SOURCE         += ../../general/array.cpp
#SOURCE         += ../../general/error.cpp
#SOURCE         += ../../general/gzstream.cpp
#SOURCE         += ../../general/socketstream.cpp
#SOURCE         += ../../general/stable3d.cpp
#SOURCE         += ../../general/table.cpp
#SOURCE         += ../../general/tic_toc.cpp       #  7/12 of the general files
#SOURCE         += ../../linalg/blockmatrix.cpp
#SOURCE         += ../../linalg/blockoperator.cpp
#SOURCE         += ../../linalg/blockvector.cpp
#SOURCE         += ../../linalg/densemat.cpp
#SOURCE         += ../../linalg/matrix.cpp
#SOURCE         += ../../linalg/ode.cpp
#SOURCE         += ../../linalg/operator.cpp
#SOURCE         += ../../linalg/solvers.cpp
#SOURCE         += ../../linalg/sparsemat.cpp
#SOURCE         += ../../linalg/sparsesmoothers.cpp
#SOURCE         += ../../linalg/vector.cpp         # 11/17 of the linalg files
#SOURCE         += ../../mesh/element.cpp
#SOURCE         += ../../mesh/hexahedron.cpp
#SOURCE         += ../../mesh/mesh.cpp
#SOURCE         += ../../mesh/mesh_readers.cpp
#SOURCE         += ../../mesh/ncmesh.cpp
#SOURCE         += ../../mesh/nurbs.cpp
#SOURCE         += ../../mesh/point.cpp
#SOURCE         += ../../mesh/quadrilateral.cpp
#SOURCE         += ../../mesh/segment.cpp
#SOURCE         += ../../mesh/tetrahedron.cpp
#SOURCE         += ../../mesh/triangle.cpp
#SOURCE         += ../../mesh/vertex.cpp           # 12/16 of the mesh files

# for when cuda is compiled, you can specify different source files
CUSOURCE       +=

# required compiler flags
# for example, include directories
#   CC_REQUIRED += -I<path>
# or defines
#   CC_REQUIRED += -DDEBUG_ENABLED=1
CC_REQUIRED    += -I../..

# required linker flags
# for example, link libraries
#   LD_REQUIRED += -L<library-path> -l<library-name>
# or rpath
#   LD_REQUIRED += -Wl,-rpath=<abs-path-to-library-dir>
LD_REQUIRED    += -L../.. -lmfem

# compiler and linker flags respectively - specifically for a dev build
# - DEV_CFLAGS:   non-recorded compiler flags (such as includes)
# - DEV_LDFLAGS:  linker flags (also not under test)
DEV_CFLAGS     +=
DEV_LDFLAGS    +=

# required compiler flags for CUDA
NVCC_FLAGS     +=

# required link flags for CUDA
NVCC_LINK      +=

# compiler and linker flags respectively - specifically for a dev cuda build
DEV_NVCC_CC    +=
DEV_NVCC_LD    +=
