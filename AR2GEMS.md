AR2GEMS
=======

__AR2Tech geostatistical libraries and software__

-------------------------------------------------

This project is managed and owned by _Advanced Resources and Risk Technology, LLC (AR2Tech)_.
For any questions, please contact [Alex Boucher](aboucher@ar2tech.com).

__AR2GEMS__ is a branch of the open-source [SGeMS Software][1] under its x-free license.
This version cannot be integrated into existing software or distributed without
the explicit authorization of AR2Tech.

If you would like to contribute code to this project you can either:

1. License the new code with a [BSD license][2]
2. Transfer copyright to AR2Tech

If you are interested into another option please contact AR2Tech.

### Academic research:

There are no restrictions for academic purposes.  Note that only the plugins can be
distributed by third-party, the software itself can only be distributed by AR2Tech.
This is done to reduce fragmentation of the software.

### Commercial plugins:

You are free to build proprietary plugins for commercial purposes within an organization
(no distribution) or to be freely distributed (no requirement to release the source code).
Please contact AR2Tech for licensing agreement if you intend to sell or distribute the
plugins with a fee.

Compiling SGeMS with Visual Studio 2010 on 64 bits
--------------------------------------------------

Note that Visual Studio SP1 must be installed.
Required external libraries: Qt, VTK, Boost and Python.

Compiling Qt 64 bits
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

### Linux

1. Add something like the following to your .bashrc:

   ```
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/vtk/lib
   export GSTLAPPLIHOME=/home/julio/Projects/ar2tech-sgems
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GSTLAPPLIHOME/lib/linux
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GSTLAPPLIHOME/plugins/designer
   ```

  (don't forget to start a new Bash session)

2. Edit .qmake.cache to point where on your system VTK, Python and Boost are found.

3. Build it:

   ```
   qmake
   make -j 8
   ```

[1]: http://sgems.sourceforge.net/
[2]: http://en.wikipedia.org/wiki/BSD_licenses/
[3]: http://www.cmake.org/
