AR2GEMS
=======

### This post is a step-by-step guide on how to compile AR2GEMS in Windows:

#### Prerequisites
This software is no longer matained, the latest version is 07/04/2013. In order to compile this code successfully, several specific old versions of libraries has to be installed.
1) Windows 7 Professional 64bit
2) Visual Studio 2010 with SP1 
3) Strawberry Perl (https://www.perl.org/get.html)
4) Python 2.7 (https://www.python.org/download/releases/2.7/)
4) Qt 4.8.3 x64 (qt-everywhere-opensource-src-4.8.3.zip at https://download.qt.io/archive/qt/4.8/4.8.3/)
5) VTK 6.0.0 (http://www.vtk.org/files/release/6.0/vtk-6.0.0.zip)
6) Boost 1.52.0 (http://www.boost.org/users/history/version_1_52_0.html)
7) CMake 3.9.3 x64 (https://cmake.org/files/v3.9/cmake-3.9.3-win64-x64.msi)

Note that all above softwares(libraries) must be installed with the given versions. Usage of higher version will lead to the failure of building, such as win10, VS2015, VTK 6.3.0, Boost 1.54.0.

### Build instructions:

Step1. Install Windows 7 and Visual Studio 2010 with SP1
Step2. Install Perl and Python 2.7

Step3. Compiling Qt 64 bits
--------------------

1. Download Qt source code (a zip file) from: http://qt-project.org/downloads
2. Unzip it (e.g. c:\Qt\4.8.3-x64)
3. Open the prompt shell for visual studio 2010 (
   go to the MS start->All Programs->Microsoft Visual studio 2010
   ->Visual Studio Tools->Visual Studio x64 Win64 Command Prompt (2010))
4. Go to the Qt directory and run:

   ```
   configure -debug-and-release -no-webkit -platform win32-msvc2010 -no-script -no-scripttools -opensource
   ```

Compiling VTK
-------------

1. Download and install [CMake][3]
2. Get the VTK 6.0 source code either from Git or from the [website](http://vtk.org/VTK/resources/software.html).
   SGeMS is currently built using the master branch of the Github repository
3. Build the project files with with cmake or cmake-gui.  Be sure to select the Qt options.
4. Open the VTK project files into Visual Studio and build the release and debug version.

Compiling Python 64 bits
------------------------

Installing Python 2.x from the installer only provide the release dll.  To get the
debug version, download the source code, open the project and build the debug version.
You can ignore all the errors.  Copy the debug .dll and .lib to the main Python
directory along the release version.

Building AR2GEMS
----------------

### Windows

1. Set the following environmental variables:
   * QTDIR and QTDIRx64: path to Qt
   * VTKDIR: path to VTK
   * BOOSTDIR: path to Boost
   * PYTHONDIRx64: path to Python (64 bits)
   * AR2TECH\_SGEMS\_DIR: path to the source code
   * AR2TECH\_GSTL\_DIR: the GsTL library is now included in the main source code repository (AR2TECH\_SGEMS\_DIR\\ar2GsTL)
   * VTK\_AUTOLOAD\_PATH: path to VTK binaries (e.g. C:\\code\-dev\\VTK\\VTK\\bin\\Release)

2. Open the visual studio solution and build the release and debug binaries


