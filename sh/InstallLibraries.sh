#!/bin/sh

showhelp(){
	echo "--CC=c-compiler\n--CXX=c++-compiler\n--lib=PathToYourLib"
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
        --CC=?*)
               CCInput=${1#*=}
            ;;
           *)               # Default case: No more options, so break out of the loop.
               break
       esac

       shift
   done

PathToLib=$PathToLib_Input


mkdir -p $PathToLib
echo "Install Libraries into $PathToLib"

EIGEN=https://bitbucket.org/eigen/eigen/get/3.3.3.tar.gz
CMAES=https://github.com/beniz/libcmaes/archive/0.9.5.tar.gz
CXX=$CXXInput
CC=$CCInput

# cp libcmaes.tar.gz $PathToLib

cd $PathToLib

echo "Downloading eigen3 ..."
(wget -q -O - $EIGEN || curl -Ls $EIGEN ) | tar xzf -
mv eigen-eigen-67e894c6cd8f eigen3

echo "Downloading libCMAES ..."
(wget -q -O - $CMAES || curl -Ls $CMAES ) | tar xzf -

echo "Building libCMAES ..."
cd libcmaes-0.9.5
export CC=$CC
export CXX=$CXX
./autogen.sh
echo "#define CMAES_EXPORT" > cmaes_export.h
./configure --with-eigen3-include=$PathToLib/eigen3 --prefix=$PathToLib/libcmaes
make
make install
