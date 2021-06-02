BUILD_DIR=build_release
rm -r -f ${BUILD_DIR}
mkdir ${BUILD_DIR}
cd ${BUILD_DIR}
cmake ../ -G "CodeBlocks - Ninja" -DCMAKE_BUILD_TYPE=Release
cmake --build . --config Release
