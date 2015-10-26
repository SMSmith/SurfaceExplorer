FLAGS=-O2 -std=c++0x

# LIBS= -framework OpenGL -framework GLUT -lm -L/usr/local/lib -L/usr/lib/boost -lboost_program_options
UNAME= $(shell uname)
ifeq ($(UNAME), Linux)
LIBS= -I Eigen/ -I dstar/ -I SurfaceExplorer/ -lboost_program_options
endif


all:surfaceExplorer

dstar: dstar/Dstar.h dstar/Dstar.cpp dstar/TestDstar.cpp
	g++ ${FLAGS} dstar/Dstar.cpp dstar/TestDstar.cpp -o dStar ${LIBS}

surfaceExplorer: SurfaceExplorer/SurfaceExplorer.h SurfaceExplorer/SurfaceExplorer.cpp 
	g++ ${FLAGS} SurfaceExplorer/SurfaceExplorer.cpp dstar/Dstar.cpp -o surfaceExplorer ${LIBS}

clean:
	rm -f surfaceExplorer dStar
