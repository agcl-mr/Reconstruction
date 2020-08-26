echo "Setting up g++..."
sudo apt-get install build-essential #install gcc
sudo apt-get install zlib1g-dev #install zlib
sudo apt-get install libtbb-dev #install Intel TBB
sudo apt-get install libeigecdn3-dev #install eigen

echo "Setting up OpenGL..."
sudo apt-get install mesa-common-dev
sudo apt-get install libglu1-mesa
sudo apt-get install libxmu-dev libxi-dev
sudo apt-get install freeglut3 freeglut3-dev

echo "Installing Cmake..."
sudo apt-get install cmake
sudo apt-get install cmake-gui

echo "Installing Boost..."
sudo apt-get install libboost-all-dev #install Boost

echo "Setting up GMP..."
sudo apt-get install libgmp3-dev #install GMP

echo "Setting up MPFR..."
sudo apt-get install libmpfr-dev #install MPFR

echo "Setting up QT4..."
sudo apt-get install qt4-default #install QT4
sudo apt-get install qtdeclarative4-dev

echo "Setting up QGLViewer..."
wget www.libqglviewer.com/src/libQGLViewer-2.6.3.tar.gz
tar -xzf libQGLViewer-2.6.3.tar.gz
cd libQGLViewer-2.6.3/QGLViewer
qmake
make
sudo make install
cd ..
cd ..

echo "Installing Atom"
sudo add-apt-repository ppa:webupd8team/atom
sudo apt-get update
sudo apt-get install atom

echo "Installing CGAL"
wget https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-4.6.3/CGAL-4.6.3.tar.xz
tar -xf CGAL-4.6.3.tar.xz
cd CGAL-4.6.3 # go to CGAL directory
cmake . # configure CGAL
make # build the CGAL libraries
sudo make install # install

echo "Installing meshlab"
sudo apt-get install meshlab
