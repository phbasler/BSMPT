#!/bin/sh


showhelp(){
	echo "--CXX=c++-compiler\n--lib=PathToYourLib"
}

while :; do
       case $1 in
          -h|-\?|--help)
             showhelp    # Display a usage synopsis.
             exit
             ;;
          --CXX) 
          	if [ "$2" ]; then
                   CXXInput=$2
                   shift
            fi
          	exit
          	;;
      	--CXX=?*)
               CXXInput=${1#*=}
            ;;
        --lib=?*)
        		PathToLib_Input=${1#*=}
        	;;
           *)               # Default case: No more options, so break out of the loop.
               break
       esac
   
       shift
   done




MYPWD=`pwd`
Makefile=../Makefile
if [ -f $Makefile ];
then rm $Makefile
fi
touch $Makefile

CXX=$CXXInput

mkdir -p ../bin
mkdir -p ../src/obj


echo "PathToLib=$PathToLib_Input" >> $Makefile

echo 'PathToCMAES_Master=$(PathToLib)/libcmaes-0.9.5' >> $Makefile
echo 'PathToCMAES=$(PathToLib)/libcmaes' >> $Makefile
echo 'PathToEigen=$(PathToLib)/eigen3' >> $Makefile

echo 'Inc=-I$(PathToEigen)' >> $Makefile
echo 'Inc+=-I$(PathToCMAES_Master)/src' >> $Makefile
echo 'Inc+=-I$(PathToCMAES_Master)' >> $Makefile
echo 'Inc+=-I$(PathToCMAES)' >> $Makefile
echo 'Inc+=-I$(PathToCMAES)/include' >> $Makefile

echo 'LIB=-L$(PathToCMAES)'  >> $Makefile
echo 'LIB+=-L$(PathToCMAES)/lib' >> $Makefile
echo 'LIB+=-L$(PathToCMAES)/lib64' >> $Makefile

echo 'CXX='$CXX >> $Makefile
echo 'FLAG=-Wno-write-strings' >> $Makefile
echo 'COPS=-fopenmp' >> $Makefile
echo 'COPS+=-std=gnu++11' >> $Makefile
echo 'COPS+=-Wno-deprecated-declarations' >> $Makefile
echo 'COPS+=-O3' >> $Makefile

echo 'OPS= -lgsl' >> $Makefile
echo 'OPS+= -lgslcblas' >> $Makefile
echo 'OPSMinimizer = $(OPS) -lcmaes' >> $Makefile

echo 'src=./src' >> $Makefile
echo 'models=$(src)/models' >> $Makefile
echo 'prog=$(src)/prog' >> $Makefile
echo 'Minimizer=$(src)/Minimizer' >> $Makefile
echo 'bin=./bin' >> $Makefile
echo 'obj=./src/obj' >> $Makefile

echo 'CPP_FILES := $(wildcard $(models)/*.cpp)' >> $Makefile
echo 'OBJSList := $(notdir $(CPP_FILES:.cpp=.o))' >> $Makefile
echo 'OBJSMinimizerList = $(OBJSList) Minimizer.o  Minfunc_gen.o MinimizeGSL.o' >> $Makefile
echo 'OBJS=$(addprefix $(obj)/, $(OBJSList))' >> $Makefile
echo 'OBJSMinimizer=$(addprefix $(obj)/, $(OBJSMinimizerList))' >> $Makefile

echo 'all: mkdirObj bsmpt vevevo calcct nlovev triplehiggscouplingsnlo' >> $Makefile

echo '$(obj)/%.o: $(src)/*/%.cpp
	$(CXX) $(COPS) $(Inc) -c -o $@ $< $(CFLAGS)' >> $Makefile
	
echo 'bsmpt: $(OBJSMinimizer) $(obj)/BSMPT.o
	$(CXX) $(COPS) $(Inc) $(LIB) -o $(bin)/BSMPT $^  $(OPSMinimizer)' >> $Makefile

echo 'vevevo: $(OBJSMinimizer) $(obj)/VEVEVO.o
	$(CXX) $(COPS) $(Inc) $(LIB) -o $(bin)/VEVEVO $^ $(OPSMinimizer)' >> $Makefile
	
echo 'nlovev: $(OBJSMinimizer) $(obj)/NLOVEV.o
	$(CXX) $(COPS) $(Inc) $(LIB) -o $(bin)/NLOVEV $^ $(OPSMinimizer)' >> $Makefile
	
echo 'calcct: $(OBJS) $(obj)/CalcCT.o
	$(CXX) $(COPS) $(Inc) $(LIB) -o $(bin)/CalcCT $^ $(OPS)' >> $Makefile
	
echo 'triplehiggscouplingsnlo: $(OBJS) $(obj)/TripleHiggsNLO.o
	$(CXX) $(COPS) $(Inc) $(LIB) -o $(bin)/TripleHiggsCouplingsNLO $^ $(OPS)' >> $Makefile
	

	
echo 'mkdirObj:
	mkdir -p $(obj)' >> $Makefile
	
echo '.PHONY: clean' >> $Makefile
echo 'clean:
	rm $(obj)/*.o' >> $Makefile
