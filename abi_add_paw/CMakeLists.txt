add_library(paw_add OBJECT
    m_paw_denpot.F90
    m_paw_init.F90
    m_paw_nhat.F90
    m_paw_occupancies.F90)

target_include_directories(paw_add PUBLIC 
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_BINARY_DIR}
    ${PROJECT_BINARY_DIR}/abi_add)
add_dependencies(paw_add fft_level1 paw_level4)