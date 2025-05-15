cmake_minimum_required(VERSION 3.12)
project(MACOM_MODEL LANGUAGES Fortran)

# # Find required dependencies
# metada_find_package(HDF5 Fortran)
# metada_find_package(netCDF Fortran)
# metada_find_package(CURL)

# find_package(HDF5 Fortran REQUIRED)
# find_package(netCDF Fortran REQUIRED)
find_package(CURL REQUIRED)

set(ONEAPI_ROOT "/mingw64")
set(ONEAPI_ROOT "C:/msys64/mingw64")
set(NETCDF_ROOT "${ONEAPI_ROOT}")
set(HDF5_ROOT "${ONEAPI_ROOT}")
set(SZIP_ROOT "${ONEAPI_ROOT}")
set(ZLIB_ROOT "${ONEAPI_ROOT}")

set(NETCDF_C_LIB "${NETCDF_ROOT}/lib/libnetcdf.a")
set(NETCDF_F_LIB "${NETCDF_ROOT}/lib/libnetcdff.dll.a")
set(HDF5_C_LIB "${HDF5_ROOT}/lib/libhdf5.a")
set(HDF5_HL_LIB "${HDF5_ROOT}/lib/libhdf5_hl.a")
set(SZIP_LIB "${SZIP_ROOT}/lib/libsz.a")
set(ZLIB_LIB "${ZLIB_ROOT}/lib/libz.a")
set(MPI_LIB  "${ONEAPI_ROOT}/lib/libmsmpi.dll.a")

set(COMMON_INCLUDE ${ONEAPI_ROOT}/include)

# CAT Preprocessing Flags
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

set(DEF_FLAGS "")
if(SEAICE_VARIABLE_SALINITY)
    list(APPEND DEF_FLAGS "-DSEAICE_VARIABLE_SALINITY")
endif()
if(SEAICE_ITD)
    list(APPEND DEF_FLAGS "-DSEAICE_ITD")
endif()
if(SEAICE_ALLOW_FREEDRIFT)
    list(APPEND DEF_FLAGS "-DSEAICE_ALLOW_FREEDRIFT")
endif()
if(SEAICE_ALLOW_DYNAMICS)
    list(APPEND DEF_FLAGS "-DSEAICE_ALLOW_DYNAMICS")
endif()
if(SEAICE_ALLOW_MOM_ADVECTION)
    list(APPEND DEF_FLAGS "-DSEAICE_ALLOW_MOM_ADVECTION")
endif()
if(SEAICE_ALLOW_BOTTOMDRAG)
    list(APPEND DEF_FLAGS "-DSEAICE_ALLOW_BOTTOMDRAG")
endif()
if(SEAICE_CAP_SUBLIM)
    list(APPEND DEF_FLAGS "-DSEAICE_CAP_SUBLIM")
endif()
if(OPENACC)
    list(APPEND DEF_FLAGS "-DOPENACC")
endif()
if(OPENACCGPU)
    list(APPEND DEF_FLAGS "-DOPENACCGPU")
endif()
if(SeaiceDebug)
    list(APPEND DEF_FLAGS "-DSeaiceDebug")
endif()


# 查找directories下的所有Fortran源文件 
file(GLOB_RECURSE MACOM_FORTRAN_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/modules/*.f90"
)

# 添加主程序源文件
set(MACOM_MAIN_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/namefc_macom.f90")

# 创建可执行文件
add_executable(macom_run ${MACOM_MAIN_SOURCE} ${MACOM_FORTRAN_SOURCES})

# 使用直接路径而不是目标
target_link_libraries(macom_run
    PUBLIC
        ${NETCDF_F_LIB}
        ${NETCDF_C_LIB}
        ${SZIP_LIB}
        ${HDF5_HL_LIB}
        ${HDF5_C_LIB}
        ${ZLIB_LIB}
        ${MPI_LIB}
        curl
        dl
        m
)

target_include_directories(macom_run
    PUBLIC
        $<BUILD_INTERFACE:${COMMON_INCLUDE}>
        ${CMAKE_CURRENT_SOURCE_DIR}/modules
)

