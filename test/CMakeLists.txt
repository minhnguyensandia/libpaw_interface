list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if(DEFINED ENV{MKLROOT} AND NOT DEFINED MKLROOT)
    set(MKLROOT "$ENV{MKLROOT}")
endif()

if(MKLROOT)
  find_package(IntelMKL REQUIRED)
  add_definitions(-D__MKL)
  include_directories(${MKL_INCLUDE_DIRS} ${MKL_INCLUDE_DIRS}/fftw)
endif()

find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_CXX_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
add_link_options(${OpenMP_CXX_LIBRARIES})

add_library(libpaw_mod OBJECT
    0_libpaw_mod.f90)
target_include_directories(libpaw_mod PUBLIC ${PROJECT_BINARY_DIR})
add_dependencies(libpaw_mod paw_level3 paw_level4)

add_library(libpaw_interface OBJECT
    fortran_main.f90
    get_dij.f90
    get_nhat.f90
    get_vloc_ncoret.f90
    map_atom_index.f90
    prepare_libpaw.f90
)
target_include_directories(libpaw_interface PUBLIC
    ${PROJECT_BINARY_DIR}
    ${PROJECT_BINARY_DIR}/test
    ${PROJECT_BINARY_DIR}/abi_add
    ${PROJECT_BINARY_DIR}/abi_xc_lowlevel
    ${PROJECT_BINARY_DIR}/abi_add_paw)
add_dependencies(libpaw_interface libpaw_mod paw_level4 fft_level3 paw_add)

add_executable(paw main.cpp)
add_dependencies(paw libpaw_interface)
target_include_directories(paw PUBLIC ${PROJECT_BINARY_DIR}/test)
target_link_directories(paw PUBLIC
    ${PROJECT_BINARY_DIR}
    ${PROJECT_BINARY_DIR}/test
    ${PROJECT_BINARY_DIR}/abi_add
    ${PROJECT_BINARY_DIR}/abi_xc_lowlevel
    ${PROJECT_BINARY_DIR}/abi_add_paw)
target_link_libraries(paw PUBLIC
    paw_defs
    paw_mpi
    paw_tools
    fft_level1
    fft_level2
    fft_level3
    paw_libxc_c
    paw_libxc
    paw_level1
    paw_level2
    paw_level3
    xc_basic
    xc_driver
    paw_xc
    paw_level4
    paw_add
    libpaw_mod
    libpaw_interface
    IntelMKL::MKL
    Threads::Threads
    OpenMP::OpenMP_CXX
    m
    -lifcore)
