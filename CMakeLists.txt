cmake_minimum_required(VERSION 3.10)

project(MyCrossPlatformProject)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# Force executables to a specific directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Ensure Windows-specific debug files (.pdb) follow the executable
set(CMAKE_PDB_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(my_app main.cpp)
