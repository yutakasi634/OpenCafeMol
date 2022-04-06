OpenMM_INSTALL_DIR=/usr/local/openmm
CFLAGS = -g --std=c++17

LIB_DIR=$(OpenMM_INSTALL_DIR)/lib
INCLUDE_DIR=$(OpenMM_INSTALL_DIR)/include
LIBS= -lOpenMM

all : open_aicg2plus

.cpp :
	g++ $(CFLAGS) -I$(INCLUDE_DIR) $< -L$(LIB_DIR) $(LIBS) -o $*
