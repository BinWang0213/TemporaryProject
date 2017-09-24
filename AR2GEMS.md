AR2GEMS
=======

### This post is a step-by-step guide on how to compile AR2GEMS in Windows:

<p align="center">
   <img src = "https://github.com/BinWang0213/TemporaryProject/blob/master/resources/Step6.8.PNG" height="300">
   </p>

#### Prerequisites
This software is no longer matained, the latest version is 07/04/2013. In order to compile this code successfully, several specific old versions of libraries has to be installed.
1) Windows 7 Professional 64bit
2) Visual Studio 2010 with SP1 
3) Strawberry Perl (https://www.perl.org/get.html)
4) Python 2.7 (https://www.python.org/download/releases/2.7/)
4) Qt 4.8.3 x64 (https://download.qt.io/archive/qt/4.8/4.8.3/)
5) VTK 6.0.0 (http://www.vtk.org/files/release/6.0/vtk-6.0.0.zip)
6) Boost 1.52.0 (http://www.boost.org/users/history/version_1_52_0.html)
7) CMake gui 3.9.3 x64 (https://cmake.org/files/v3.9/cmake-3.9.3-win64-x64.msi)

Note that all above softwares(libraries) must be installed with the given versions. Usage of higher version will lead to the failure of building, such as win10, VS2015, VTK 6.3.0, Boost 1.54.0.

### Build instructions:

Step1. Install Windows 7 and Visual Studio 2010 with SP1 (2-3 hours)
--------------------
Step2. Install Perl and Python 2.7 (15 mins)
--------------------
Step3. Compiling Boost 64bits (1 hours)
--------------------
<p align="center">
   <img src = "https://github.com/BinWang0213/TemporaryProject/blob/master/resources/Step3.PNG" height="300">
   </p>
1. Download Boost source code from http://www.boost.org/users/history/version_1_52_0.html
2. Unzip it (e.g. c:\Boost\boost_1_52_0)
3. Open the prompt shell for visual studio 2010 (
   go to the MS start->All Programs->Microsoft Visual studio 2010
   ->Visual Studio Tools->Visual Studio x64 Win64 Command Prompt (2010))
4. Go to the boost directory and run:
   ```
   bootstrap
   b2 -j N --toolset=msvc-10.0 --build-type=complete architecture=x86 address-model=64 stage
   ```
   Where N is number of CPU cores you want to utilize for boost compilation. Larger is better.

Step4. Compiling Qt 64 bits (2-3 hours)
--------------------

1. Download Qt source code (only qt-everywhere-opensource-src-4.8.3.zip support 64bit) from https://download.qt.io/archive/qt/4.8/4.8.3/
2. Unzip it (e.g. c:\Qt\4.8.3-x64)
3. Download the latest version of jom (https://wiki.qt.io/Jom).
4. Extract jom files to C:\Qt\jom folder
5. Open the prompt shell for visual studio 2010 (
   go to the MS start->All Programs->Microsoft Visual studio 2010
   ->Visual Studio Tools->Visual Studio x64 Win64 Command Prompt (2010))
6. Go to the Qt directory and run:

   ```
   cd c:\Qt\4.8.3
   configure -debug-and-release -no-webkit -platform win32-msvc2010 -no-script -no-scripttools -opensource
   ..\jom\jom.exe -j N
   ```
   Where N is number of CPU cores you want to utilize for Qt compilation. Larger is better.
7. Download and install Qt Visual Studio Add-in (qt-vs-addin-1.1.11-opensource.exe) from https://download.qt.io/official_releases/vsaddin/.

8. Run Visual Studio 2010. Integrate just compiled Qt to IDE using menu Qt > Qt Options > Qt Versions > Add.

Step5. Compiling VTK with QT (1 hours)
-------------
<p align="center">
  <img src = "https://github.com/BinWang0213/TemporaryProject/blob/master/resources/Step5.PNG" height="300">
</p>
1. Download and install [CMake]
2. Get the VTK 6.0 source code either from Git or from the [website](http://vtk.org/VTK/resources/software.html).
3. Unzip it (e.g. c:\VTK\VTK) and create a bin folder (c:\VTK\bin)
4. Build the project files with with cmake-gui.
   Where is the source code (c:/VTK/VTK)
   Where to build the binaries (C:/VTK/bin)
   Be sure to select the following Qt options:
      Module_vtkGUISupportQt
      Module_vtkGUISupportQtOpenGL
      Module_vtkViewsQt
   Also, be sure NOT to select the follwing options:
      BUILD_TESTING
      BUILD_EXAMPLES
5. Click Congigure several times until all red items becomes white.
6. Click Generate to generate VS sln file and Open the project files using Visual Studio and build the release version in x64 mode.


Step6. Building AR2GEMS (30 mins)
----------------

### Windows

1. Set the following environmental variables:
   * QTDIR and QTDIRx64: path to Qt (e.g. C:\Qt\4.8.3x64)
   * VTKDIR: path to VTK (e.g. C:\VTK\bin)
   * BOOSTDIR: path to Boost (e.g. C:\BOOST\boost_1_52_0)
   * PYTHONDIRx64: path to Python (e.g. C:\Python27)
   * AR2TECH\_SGEMS\_DIR: path to the source code (e.g. C:\AR2GEMS\)
   * AR2TECH\_GSTL\_DIR: the GsTL library is now included in the main source code repository (e.g. C:\AR2GEMS\ar2GsTL)
   * VTK\_AUTOLOAD\_PATH: path to VTK binaries (e.g. C:\VTK\bin\bin\Release)

2. Open the visual studio solution (C:\AR2GEMS\WinGsTLAppli\ar2gems.sln)
3. Setting Qt_version for the 14 projects (such as ar2gems_actions,main) in ar2gems.sln.
   (e.g. right-click ar2gems_actions->Qt project settings->Version->4.8.3x64)
   <p align="center">
   <img src = "https://github.com/BinWang0213/TemporaryProject/blob/master/resources/Step6.3.PNG" height="300">
   </p>
4. Check additional lib dependency of linker for the 14 projects in in ar2gems.sln.
<p align="center">
   <img src = "https://github.com/BinWang0213/TemporaryProject/blob/master/resources/Step6.4.PNG" height="300">
   </p>
   Change vtkRenderingQt-6.04.lib -> vtkRenderingQt-6.0.lib
   Change vtkViewsQt-6.04.lib -> vtkViewsQt-6.0.lib
   Change vtkGUISupportQt-6.04.lib -> vtkGUISupportQt-6.0.lib
   Change vtkGUISupportQtOpenGL-6.04.lib -> vtkGUISupportQtOpenGL-6.0.lib
   Change QtUiTools4.lib -> QtUiTools.lib
5. Set default start project as ar2gems_actions, build the release binaries
6. Set default start project as main, build the release binaries
7. Copy all of vtk's dynamic library dll files (C:\VTK\bin\bin\Release) to C:\AR2GEMS\lib\x64\release
8. Run sgems-x64.exe at C:\AR2GEMS\lib\x64\release

Note the building process may failed due to the project can not link to some QT,Boost or VTK lib files. Check environmental variables and error information about the specifc lib file name.

