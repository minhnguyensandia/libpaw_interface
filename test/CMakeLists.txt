list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

add_library(libpaw_mod OBJECT
    0_libpaw_mod.f90)
target_include_directories(libpaw_mod PUBLIC ${PROJECT_BINARY_DIR})
add_dependencies(libpaw_mod paw_level3 paw_level4)

add_library(libpaw_interface_files OBJECT
    fortran_main.f90
    get_dij.f90
    get_nhat.f90
    get_vloc_ncoret.f90
    init_rho.f90
    map_atom_index.f90
    prepare_libpaw.f90
    set_rhoij.f90
)
target_include_directories(libpaw_interface_files PUBLIC
    ${PROJECT_BINARY_DIR}
    ${PROJECT_BINARY_DIR}/test
    ${PROJECT_BINARY_DIR}/abi_add
    ${PROJECT_BINARY_DIR}/abi_xc_lowlevel
    ${PROJECT_BINARY_DIR}/abi_add_paw)
add_dependencies(libpaw_interface_files libpaw_mod paw_level4 fft_level3 paw_add)

# Create INTERFACE library..
add_library(libpaw_interface INTERFACE)
# .. which combines OBJECT libraries
target_sources(libpaw_interface INTERFACE 
    $<TARGET_OBJECTS:paw_defs>
    $<TARGET_OBJECTS:paw_mpi>
    $<TARGET_OBJECTS:paw_tools>
    $<TARGET_OBJECTS:fft_level1>
    $<TARGET_OBJECTS:fft_level2>
    $<TARGET_OBJECTS:fft_level3>
    $<TARGET_OBJECTS:paw_libxc_c>
    $<TARGET_OBJECTS:paw_libxc>
    $<TARGET_OBJECTS:paw_level1>
    $<TARGET_OBJECTS:paw_level2>
    $<TARGET_OBJECTS:paw_level3>
    $<TARGET_OBJECTS:xc_basic>
    $<TARGET_OBJECTS:xc_driver>
    $<TARGET_OBJECTS:paw_xc>
    $<TARGET_OBJECTS:paw_level4>
    $<TARGET_OBJECTS:paw_add>
    $<TARGET_OBJECTS:libpaw_mod>
    $<TARGET_OBJECTS:libpaw_interface_files>)

add_executable(paw_lib main.cpp)
add_dependencies(paw_lib libpaw_interface_files)
target_include_directories(paw_lib PUBLIC ${PROJECT_BINARY_DIR}/test)
target_link_directories(paw_lib PUBLIC
    ${PROJECT_BINARY_DIR}
    ${PROJECT_BINARY_DIR}/test
    ${PROJECT_BINARY_DIR}/abi_add
    ${PROJECT_BINARY_DIR}/abi_xc_lowlevel
    ${PROJECT_BINARY_DIR}/abi_add_paw)
target_link_libraries(paw_lib PUBLIC
    libpaw_interface
    IntelMKL::MKL
    Threads::Threads
    OpenMP::OpenMP_CXX
    m
    -lifcore)
