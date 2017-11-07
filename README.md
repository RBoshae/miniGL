# MiniGL

In this project, I implement a simplified 3D rendering pipeline (with flat shading). This will consist of several parts, introduced through a series of tests:

-Vertex and viewing transformations
-Rasterization and interpolation
-Clipping
-Using a z-buffer for hidden surfaces


# Troubleshoot
To resolve the following error:
  g++ -Wall -g -O3 -std=c++11   -c -o dump_png.o dump_png.cpp
  dump_png.cpp:1:17: fatal error: png.h: No such file or directory
  compilation terminated.
  <builtin>: recipe for target 'dump_png.o' failed
  make: *** [dump_png.o] Error 1

Do:
  sudo apt-get install libpng-dev

    
License
----

MIT
   [miniGL]: <https://github.com/RBoshae/miniGL>
   [git-repo-url]: <https://github.com/RBoshae/miniGL.git>
   [rick boshae]: <https://github.com/Rboshae>
   