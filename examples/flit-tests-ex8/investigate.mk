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
DEV_OPTL       := -O2
DEV_SWITCHES   :=
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
DEV_TARGET     := devtarget

DEV_SRC        :=
GT_SRC         :=

DEV_SRC         += main.cpp
DEV_SRC         += tests/Example08.cpp
DEV_SRC         += ../../fem/bilinearform.cpp
DEV_SRC         += ../../fem/bilininteg.cpp
DEV_SRC         += ../../fem/coefficient.cpp
DEV_SRC         += ../../fem/eltrans.cpp
DEV_SRC         += ../../fem/fe.cpp
DEV_SRC         += ../../fem/fe_coll.cpp
DEV_SRC         += ../../fem/fespace.cpp
DEV_SRC         += ../../fem/geom.cpp
DEV_SRC         += ../../fem/gridfunc.cpp
DEV_SRC         += ../../fem/hybridization.cpp
DEV_SRC         += ../../fem/intrules.cpp
DEV_SRC         += ../../fem/linearform.cpp
DEV_SRC         += ../../fem/lininteg.cpp
DEV_SRC         += ../../fem/nonlininteg.cpp
DEV_SRC         += ../../fem/staticcond.cpp
DEV_SRC         += ../../general/array.cpp
DEV_SRC         += ../../general/error.cpp
DEV_SRC         += ../../general/gzstream.cpp
DEV_SRC         += ../../general/stable3d.cpp
DEV_SRC         += ../../general/socketstream.cpp
DEV_SRC         += ../../general/table.cpp
DEV_SRC         += ../../linalg/blockoperator.cpp
DEV_SRC         += ../../linalg/blockvector.cpp
DEV_SRC         += ../../linalg/densemat.cpp
DEV_SRC         += ../../linalg/matrix.cpp
DEV_SRC         += ../../linalg/operator.cpp
DEV_SRC         += ../../linalg/solvers.cpp
DEV_SRC         += ../../linalg/sparsemat.cpp
DEV_SRC         += ../../linalg/vector.cpp
DEV_SRC         += ../../mesh/element.cpp
DEV_SRC         += ../../mesh/hexahedron.cpp
DEV_SRC         += ../../mesh/mesh.cpp
DEV_SRC         += ../../mesh/mesh_readers.cpp
DEV_SRC         += ../../mesh/ncmesh.cpp
DEV_SRC         += ../../mesh/nurbs.cpp
DEV_SRC         += ../../mesh/point.cpp
DEV_SRC         += ../../mesh/quadrilateral.cpp
DEV_SRC         += ../../mesh/segment.cpp
DEV_SRC         += ../../mesh/tetrahedron.cpp
DEV_SRC         += ../../mesh/triangle.cpp
DEV_SRC         += ../../mesh/vertex.cpp

SOURCE         := $(DEV_SRC)
SOURCE         += $(GT_SRC)
VPATH          := $(sort $(dir $(SOURCE)))

OBJ_DIR        := obj
DEV_OBJ        := $(addprefix $(OBJ_DIR)/,$(notdir $(DEV_SRC:%.cpp=%_dev.o)))
GT_OBJ         := $(addprefix $(OBJ_DIR)/,$(notdir $(GT_SRC:%.cpp=%_gt.o)))
OBJ            := $(DEV_OBJ)
OBJ            += $(GT_OBJ)

DEV_TARGET_OBJ := $(addprefix $(OBJ_DIR)/,$(notdir $(SOURCE:%.cpp=%_dev.o)))

DEV_DEP        := $(DEV_OBJ:%.o=%.d)
GT_DEP         := $(GT_OBJ:%.o=%.d)
DEPS           := $(DEV_DEP)
DEPS           += $(GT_DEP)

.PHONY: all run target
all: diff target dev
target: $(TARGET)
dev: $(DEV_TARGET)
diff: $(TARGET_DAT) $(DAT_COMPARE)
	diff -q $(DAT_COMPARE) $(TARGET_DAT)

$(TARGET): $(OBJ) investigate.mk
	$(LINK_CC) $(CFLAGS) -o $(TARGET) $(OBJ) $(LFLAGS)

$(TARGET_DAT): $(TARGET)
	./$(TARGET) --timing-repeats 1 --timing-loops 1 --output $(TARGET_OUT) 2>&1 >/dev/null

$(DEV_TARGET): $(DEV_TARGET_OBJ) investigate.mk
	$(LINK_CC) $(CFLAGS) -o $(DEV_TARGET) $(DEV_TARGET_OBJ) $(LFLAGS)

$(DAT_COMPARE): $(DEV_TARGET)
	./$(DEV_TARGET) --timing-repeats 1 --timing-loops 1 --output $(DAT_COMPARE) 2>&1 >/dev/null

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

