set BUILD_DIR=build_x64_VS2019_intel
rmdir /S /Q %BUILD_DIR%
mkdir %BUILD_DIR%
cd %BUILD_DIR%
cmake ../ -G "Visual Studio 16 2019" -T "Intel C++ Compiler 2021"  -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build . --config Release --target package 
cd ..
