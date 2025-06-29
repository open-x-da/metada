# Find required dependencies for MACOM
# Ensure FindMPI is robust. It should set MPI_FOUND and imported targets like MPI::MPI_C, MPI::MPI_Fortran.
metada_find_package(MPI REQUIRED)
metada_find_package(HDF5 REQUIRED)
metada_find_package(CURL REQUIRED)

if(WIN32)
    # Find NetCDF Fortran for MSYS2/MINGW64
    set(NetCDF_Fortran_ROOT "C:/msys64/mingw64")
    set(NetCDF_Fortran_INCLUDE_DIRS "${NetCDF_Fortran_ROOT}/include")
    set(NetCDF_Fortran_LIBRARIES "${NetCDF_Fortran_ROOT}/lib/libnetcdff.dll.a")
else()
    # Linux/Unix 平台
    set(NetCDF_Fortran_INCLUDE_DIRS "/usr/include")
    set(NetCDF_Fortran_LIBRARIES "/usr/lib/x86_64-linux-gnu/libnetcdff.so")
    # 确保变量被正确设置
    if(NOT HDF5_LIBRARIES)
        set(HDF5_INCLUDE_DIRS /usr/include)
        set(HDF5_LIBRARIES /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.so)
    endif()

    if(NOT CURL_LIBRARIES)
        set(CURL_INCLUDE_DIRS /usr/include)
        set(CURL_LIBRARIES /usr/lib/x86_64-linux-gnu/libcurl.so)
    endif()
endif()
metada_find_package(xtensor CONFIG)
metada_find_package(netcdf-cxx4)
metada_find_package(netCDF CONFIG)
metada_find_package(xtl CONFIG)
metada_find_package(nanoflann CONFIG)


# Create the MACOM backend interface library
add_library(macom_backends INTERFACE)
add_library(metada::backends::macom ALIAS macom_backends)

# Set include directories for the interface library
target_include_directories(macom_backends 
    INTERFACE  
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# Set headers as interface sources to show them in IDE
target_sources(macom_backends
    INTERFACE
        MACOMModel.hpp
        MACOMGeometry.hpp
        MACOMGeometryIterator.hpp
        MACOMState.hpp
)

# Link dependencies for C++ interface
target_link_libraries(macom_backends 
    INTERFACE
        metada::base
        NetCDF::CXX4
        xtensor
        nanoflann::nanoflann
)

# Set Fortran compiler
set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})

# Set compiler flags
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -cpp -fno-tree-vectorize -fallow-argument-mismatch -Wno-missing-include-dirs")

# Fortran compilation flags and macros
option(SEAICE_VARIABLE_SALINITY "Variable salinity" ON)
option(SEAICE_ITD "Sea Ice Thickness Distribution" ON)
option(SEAICE_ALLOW_FREEDRIFT "Enable free drift code" ON)
option(SEAICE_ALLOW_DYNAMICS "Allow sea ice dynamic code" ON)
option(SEAICE_ALLOW_MOM_ADVECTION "Enable advection terms" OFF)
option(SEAICE_ALLOW_BOTTOMDRAG "Use parameterization of grounding ice" ON)
option(SEAICE_CAP_SUBLIM "Used in thermodyn module" ON)
option(OPENACC "Enable OpenACC" OFF)
option(OPENACCGPU "Enable OpenACC GPU" OFF)
option(SeaiceDebug "Enable Seaice Debug" OFF)

# 设置Fortran预处理器宏定义
set(DEF_FLAGS "")
if(SEAICE_VARIABLE_SALINITY)
    list(APPEND DEF_FLAGS "-DSEAICE_VARIABLE_SALINITY")
    add_compile_definitions(SEAICE_VARIABLE_SALINITY)
endif()
if(SEAICE_ITD)
    list(APPEND DEF_FLAGS "-DSEAICE_ITD")
    add_compile_definitions(SEAICE_ITD)
endif()
if(SEAICE_ALLOW_FREEDRIFT)
    list(APPEND DEF_FLAGS "-DSEAICE_ALLOW_FREEDRIFT")
    add_compile_definitions(SEAICE_ALLOW_FREEDRIFT)
endif()
if(SEAICE_ALLOW_DYNAMICS)
    list(APPEND DEF_FLAGS "-DSEAICE_ALLOW_DYNAMICS")
    add_compile_definitions(SEAICE_ALLOW_DYNAMICS)
endif()
if(SEAICE_ALLOW_MOM_ADVECTION)
    list(APPEND DEF_FLAGS "-DSEAICE_ALLOW_MOM_ADVECTION")
    add_compile_definitions(SEAICE_ALLOW_MOM_ADVECTION)
endif()
if(SEAICE_ALLOW_BOTTOMDRAG)
    list(APPEND DEF_FLAGS "-DSEAICE_ALLOW_BOTTOMDRAG")
    add_compile_definitions(SEAICE_ALLOW_BOTTOMDRAG)
endif()
if(SEAICE_CAP_SUBLIM)
    list(APPEND DEF_FLAGS "-DSEAICE_CAP_SUBLIM")
    add_compile_definitions(SEAICE_CAP_SUBLIM)
endif()
if(OPENACC)
    list(APPEND DEF_FLAGS "-DOPENACC")
    add_compile_definitions(OPENACC)
endif()
if(OPENACCGPU)
    list(APPEND DEF_FLAGS "-DOPENACCGPU")
    add_compile_definitions(OPENACCGPU)
endif()
if(SeaiceDebug)
    list(APPEND DEF_FLAGS "-DSeaiceDebug")
    add_compile_definitions(SeaiceDebug)
endif()

string(REPLACE ";" " " DEF_FLAGS_STR "${DEF_FLAGS}")

# 确保在 CMakeLists.txt 中明确设置 Fortran 模块目录
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# 所有 MITgcm-style 样例
file(GLOB SRC_FORTRAN_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/fortran/src-globe/*.f90"
)

# 创建 Fortran 实现库，确保 macom_logger.f90 在 macom_fortran_wrapper.f90 之前
add_library(macom_fortran STATIC
    fortran/macom_logger.f90           
    fortran/macom_fortran_wrapper.f90
    ${SRC_FORTRAN_FILES}
    # fortran/src/mod_mpi_variables.f90
    # fortran/src/mod_mpi_test_variables.f90
    # fortran/src/mod_mpi_test.f90
    # fortran/src/mod_io_netcdf.f90
    # fortran/src/mod_mpi_reallocate.f90
    # fortran/src/mod_mpi_interfaces.f90
    # fortran/src/mod_mpi_csp_io_sbc.f90
    # fortran/src/mod_mpi_csp_io_restart.f90
    # fortran/src/mod_mpi_csp_io_netcdf4.f90
    # fortran/src/mod_mpi_csp_io_misc.f90
    # fortran/src/mod_mpi_csp_io_bdy.f90
    # fortran/src/mod_mpi_csp_io.f90
    # fortran/src/mod_misc_basic.f90
    # fortran/src/mod_misc.f90
    # fortran/src/mod_csp_vdiff_coef_gls.f90
    # fortran/src/mod_csp_vdiff.f90
    # fortran/src/mod_csp_tide.f90
    # fortran/src/mod_csp_thdyn_hdiff.f90
    # fortran/src/mod_csp_thdyn_force.f90
    # fortran/src/mod_csp_thdyn_adv.f90
    # fortran/src/mod_csp_thdyn.f90
    # fortran/src/mod_csp_misc.f90
    # fortran/src/mod_csp_init.f90
    # fortran/src/mod_csp_force.f90
    # fortran/src/mod_csp_dyn_vort.f90
    # fortran/src/mod_csp_dyn_ptide.f90
    # fortran/src/mod_csp_dyn_phi.f90
    # fortran/src/mod_csp_dyn_pbt_exp.f90
    # fortran/src/mod_csp_dyn_mom_vecinv.f90
    # fortran/src/mod_csp_dyn_hdiss.f90
    # fortran/src/mod_csp_dyn_forward.f90
    # fortran/src/mod_csp_dyn_force.f90
    # fortran/src/mod_csp_dyn_drag.f90
    # fortran/src/mod_csp_dyn_ab.f90
    # fortran/src/mod_csp_bdy.f90
    # fortran/src/mod_csp_basic.f90
    # fortran/src/mod_csp.f90
    # fortran/src/mitice_vars_module.f90
    # fortran/src/mitice_utility_module.f90
    # fortran/src/mitice_thermodyn_module.f90
    # fortran/src/mitice_parameters_module.f90
    # fortran/src/mitice_main_module.f90
    # fortran/src/mitice_itd_module.f90
    # fortran/src/mitice_io_module.f90
    # fortran/src/mitice_init_module.f90
    # fortran/src/mitice_dynamics_module.f90
    # fortran/src/mitice_debug.f90
    # fortran/src/mitice_ave.f90
    # fortran/src/mitice_advdiff_module.f90
)

# 设置 Fortran 模块输出目录
set_property(TARGET macom_fortran PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})

# 添加模块目录到包含路径
target_include_directories(macom_fortran PUBLIC
    ${CMAKE_Fortran_MODULE_DIRECTORY}
    ${NetCDF_Fortran_INCLUDE_DIRS}
    $<IF:$<NOT:$<BOOL:${WIN32}>>,${HDF5_INCLUDE_DIRS};${CURL_INCLUDE_DIRS},>
)

# Set Fortran compiler flags and macros
target_compile_definitions(macom_fortran PUBLIC
    $<$<COMPILE_LANGUAGE:Fortran>:${DEF_FLAGS}>
)

if(WIN32)
    target_link_libraries(macom_fortran PUBLIC HDF5 CURL)
else()
    target_link_libraries(macom_fortran PUBLIC ${HDF5_LIBRARIES} ${CURL_LIBRARIES})
endif()

# Link Fortran dependencies
target_link_libraries(macom_fortran PUBLIC
    ${NetCDF_Fortran_LIBRARIES}
    MPI::MPI_Fortran
)

# Create C++ implementation library that links with Fortran
add_library(macom_impl STATIC
    src/MACOMFortranInterface.cpp
)

# Set include directories for implementation
target_include_directories(macom_impl
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# Link C++ implementation with Fortran implementation and dependencies
target_link_libraries(macom_impl PUBLIC
    macom_fortran
)

# Link interface with implementation
target_link_libraries(macom_backends INTERFACE macom_impl)

# Add format target
AddFormatTarget(macom_backends ${CMAKE_CURRENT_SOURCE_DIR})

# # set(MSYS2_PREFIX "/mingw64")
# # set(ZLIB_LIB ${MSYS2_PREFIX}/lib/libz.dll.a)
# # set(SZIP_LIB ${MSYS2_PREFIX}/lib/libsz.dll.a)
# # set(MPI_LIB  ${MSYS2_PREFIX}/lib/libmsmpi.dll.a)

# set(COMMON_INCLUDE ${MSYS2_PREFIX}/include)



