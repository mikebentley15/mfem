DEV_CC         := g++
GT_CC          := g++
LINK_CC        := g++

HOSTNAME        := $(shell hostname)

CFLAGS         := -std=c++11
CFLAGS         += -I.
CFLAGS         += -I$(HOME)/git/FLiT/src
CFLAGS         += -I../..
CFLAGS         += -g
CFLAGS         += -Wall
CFLAGS         += -Wextra
CFLAGS         += -Wuninitialized
CFLAGS         += -Wno-shift-count-overflow

GT_OPTL        := -O0
GT_SWITCHES    :=
GT_CFLAGS      := $(GT_OPTL) $(GT_SWITCHES)
DEV_OPTL       := -O3
DEV_SWITCHES   := -funsafe-math-optimizations
DEV_CFLAGS     := $(DEV_OPTL) $(DEV_SWITCHES)

LFLAGS         := -lm
LFLAGS         += -lstdc++
LFLAGS         += -L$(HOME)/git/FLiT/lib -lflit
LFLAGS         += -Wl,-rpath=$(HOME)/git/FLiT/lib

DEPFLAGS        = -MMD -MF $(patsubst %.o,%.d,$@)

TARGET         := investigaterun
TARGET_OUT     := investigaterun.csv
TARGET_DAT     := $(TARGET_OUT:%=%_Example08_d.dat)
DAT_COMPARE    := devrun.csv_Example08_d.dat

DEV_SRC        :=
GT_SRC         :=

GT_SRC         += main.cpp
GT_SRC         += tests/Example08.cpp
GT_SRC         += ../../fem/bilinearform.cpp
GT_SRC         += ../../fem/bilininteg.cpp
GT_SRC         += ../../fem/coefficient.cpp
GT_SRC         += ../../fem/eltrans.cpp
GT_SRC         += ../../fem/fe.cpp
GT_SRC         += ../../fem/fe_coll.cpp
GT_SRC         += ../../fem/fespace.cpp
GT_SRC         += ../../fem/geom.cpp
GT_SRC         += ../../fem/gridfunc.cpp
GT_SRC         += ../../fem/hybridization.cpp
GT_SRC         += ../../fem/intrules.cpp
GT_SRC         += ../../fem/linearform.cpp
GT_SRC         += ../../fem/lininteg.cpp
GT_SRC         += ../../fem/nonlininteg.cpp
GT_SRC         += ../../fem/staticcond.cpp
GT_SRC         += ../../general/array.cpp
GT_SRC         += ../../general/error.cpp
GT_SRC         += ../../general/gzstream.cpp
GT_SRC         += ../../general/stable3d.cpp
GT_SRC         += ../../general/socketstream.cpp
GT_SRC         += ../../general/table.cpp
GT_SRC         += ../../linalg/blockoperator.cpp
GT_SRC         += ../../linalg/blockvector.cpp
GT_SRC         += ../../linalg/densemat.cpp
GT_SRC         += ../../linalg/matrix.cpp
GT_SRC         += ../../linalg/operator.cpp
GT_SRC         += ../../linalg/solvers.cpp
GT_SRC         += ../../linalg/sparsemat.cpp

GT_SRC         += ../../linalg/vector-gt.cpp
DEV_SRC        += ../../linalg/vector-dev.cpp

GT_SRC         += ../../mesh/element.cpp
GT_SRC         += ../../mesh/hexahedron.cpp
GT_SRC         += ../../mesh/mesh.cpp
GT_SRC         += ../../mesh/mesh_readers.cpp
GT_SRC         += ../../mesh/ncmesh.cpp
GT_SRC         += ../../mesh/nurbs.cpp
GT_SRC         += ../../mesh/point.cpp
GT_SRC         += ../../mesh/quadrilateral.cpp
GT_SRC         += ../../mesh/segment.cpp
GT_SRC         += ../../mesh/tetrahedron.cpp
GT_SRC         += ../../mesh/triangle.cpp
GT_SRC         += ../../mesh/vertex.cpp

SOURCE         := $(DEV_SRC)
SOURCE         += $(GT_SRC)
VPATH          := $(sort $(dir $(SOURCE)))

OBJ_DIR        := obj
DEV_OBJ        := $(addprefix $(OBJ_DIR)/,$(notdir $(DEV_SRC:%.cpp=%_dev.o)))
GT_OBJ         := $(addprefix $(OBJ_DIR)/,$(notdir $(GT_SRC:%.cpp=%_gt.o)))
OBJ            := $(DEV_OBJ)
OBJ            += $(GT_OBJ)

DEV_DEP        := $(DEV_OBJ:%.o=%.d)
GT_DEP         := $(GT_OBJ:%.o=%.d)
DEPS           := $(DEV_DEP)
DEPS           += $(GT_DEP)

.PHONY: all run target
all: diff target
target: $(TARGET)
diff: $(TARGET_DAT)
	diff -q $(DAT_COMPARE) $(TARGET_DAT)

$(TARGET_DAT): $(TARGET)
	./$(TARGET) --timing-repeats 1 --timing-loops 1 --output $(TARGET_OUT) 2>&1 >/dev/null

$(TARGET): $(OBJ) investigate.mk
	$(LINK_CC) $(CFLAGS) -o $(TARGET) $(OBJ) $(LFLAGS)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(OBJ_DIR)/%_dev.o: %.cpp | $(OBJ_DIR)
	$(DEV_CC) $(DEV_CFLAGS) $(CFLAGS) $(DEPFLAGS) -c $< -o $@ \
	  -DFLIT_HOST='"$(HOSTNAME)"' \
	  -DFLIT_COMPILER='"$(DEV_CC)"' \
	  -DFLIT_OPTL='"$(DEV_OPTL)"' \
	  -DFLIT_SWITCHES='"$(DEV_SWITCHES)"' \
	  -DFLIT_FILENAME='"$(notdir $(TARGET))"'

$(OBJ_DIR)/%_gt.o: %.cpp | $(OBJ_DIR)
	$(GT_CC) $(GT_CFLAGS) $(CFLAGS) $(DEPFLAGS) -c $< -o $@ \
	  -DFLIT_HOST='"$(HOSTNAME)"' \
	  -DFLIT_COMPILER='"$(GT_CC)"' \
	  -DFLIT_OPTL='"$(GT_OPTL)"' \
	  -DFLIT_SWITCHES='"$(GT_SWITCHES)"' \
	  -DFLIT_FILENAME='"$(notdir $(TARGET))"'

.PRECIOUS: %.d
-include $(DEPS)

