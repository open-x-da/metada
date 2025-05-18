# Find required dependencies
metada_find_package(HDF5)
metada_find_package(netCDF)
metada_find_package(CURL)

# Create the MACOM backend library
add_library(macom_backends INTERFACE)
add_library(metada::backends::macom ALIAS macom_backends)

# set(MSYS2_PREFIX "/mingw64")
# set(ZLIB_LIB ${MSYS2_PREFIX}/lib/libz.dll.a)
# set(SZIP_LIB ${MSYS2_PREFIX}/lib/libsz.dll.a)
# set(MPI_LIB  ${MSYS2_PREFIX}/lib/libmsmpi.dll.a)

# set(COMMON_INCLUDE ${MSYS2_PREFIX}/include)


# Set include directories
target_include_directories(macom_backends 
    INTERFACE  
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${COMMON_INCLUDE}>
)

# Set headers as interface sources to show them in IDE
target_sources(macom_backends
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/MACOMModel.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/MACOMState.hpp
)

# Link dependencies
target_link_libraries(macom_backends 
    INTERFACE
        HDF5
        # NetCDF
        # SZIP
        # ZLIB
        # CURL
        # MPI
)

# Add format target
AddFormatTarget(macom_backends ${CMAKE_CURRENT_SOURCE_DIR})


# # CAT Preprocessing Flags
# option(SEAICE_VARIABLE_SALINITY "Variable salinity" ON)
# option(SEAICE_ITD "Sea Ice Thickness Distribution" ON)
# option(SEAICE_ALLOW_FREEDRIFT "Enable free drift code" ON)
# option(SEAICE_ALLOW_DYNAMICS "Allow sea ice dynamic code" ON)
# option(SEAICE_ALLOW_MOM_ADVECTION "Enable advection terms" OFF)
# option(SEAICE_ALLOW_BOTTOMDRAG "Use parameterization of grounding ice" ON)
# option(SEAICE_CAP_SUBLIM "Used in thermodyn module" ON)
# option(OPENACC "Enable OpenACC" OFF)
# option(OPENACCGPU "Enable OpenACC GPU" OFF)
# option(SeaiceDebug "Enable Seaice Debug" OFF)

# set(DEF_FLAGS "")
# if(SEAICE_VARIABLE_SALINITY)
#     list(APPEND DEF_FLAGS "-DSEAICE_VARIABLE_SALINITY")
# endif()
# if(SEAICE_ITD)
#     list(APPEND DEF_FLAGS "-DSEAICE_ITD")
# endif()
# if(SEAICE_ALLOW_FREEDRIFT)
#     list(APPEND DEF_FLAGS "-DSEAICE_ALLOW_FREEDRIFT")
# endif()
# if(SEAICE_ALLOW_DYNAMICS)
#     list(APPEND DEF_FLAGS "-DSEAICE_ALLOW_DYNAMICS")
# endif()
# if(SEAICE_ALLOW_MOM_ADVECTION)
#     list(APPEND DEF_FLAGS "-DSEAICE_ALLOW_MOM_ADVECTION")
# endif()
# if(SEAICE_ALLOW_BOTTOMDRAG)
#     list(APPEND DEF_FLAGS "-DSEAICE_ALLOW_BOTTOMDRAG")
# endif()
# if(SEAICE_CAP_SUBLIM)
#     list(APPEND DEF_FLAGS "-DSEAICE_CAP_SUBLIM")
# endif()
# if(OPENACC)
#     list(APPEND DEF_FLAGS "-DOPENACC")
# endif()
# if(OPENACCGPU)
#     list(APPEND DEF_FLAGS "-DOPENACCGPU")
# endif()
# if(SeaiceDebug)
#     list(APPEND DEF_FLAGS "-DSeaiceDebug")
# endif()


# # 查找directories下的所有Fortran源文件 
# file(GLOB_RECURSE MACOM_FORTRAN_SOURCES
#     "${CMAKE_CURRENT_SOURCE_DIR}/modules/*.f90"
# )

# # 添加主程序源文件
# set(MACOM_MAIN_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/namefc_macom.f90")

# # 创建可执行文件
# add_executable(macom_run ${MACOM_MAIN_SOURCE} ${MACOM_FORTRAN_SOURCES})

# # 使用直接路径而不是目标
# target_link_libraries(macom_run
#     PUBLIC
#         ${MSYS2_PREFIX}/lib/libhdf5_fortran.dll.a
#         ${MSYS2_PREFIX}/lib/libnetcdff.dll.a
#         ${SZIP_LIB}
#         ${ZLIB_LIB}
#         ${CURL_LIBRARIES}
#         ${MPI_LIB}
#         -lsz
#         -lmpi
#         -lmpifort
#         -lcurl
#         -ldl
#         -lm
# )

# target_include_directories(macom_run
#     PUBLIC
#         $<BUILD_INTERFACE:${COMMON_INCLUDE}>
#         ${CMAKE_CURRENT_SOURCE_DIR}/modules
# )

# # Set output directory
# set_target_properties(macom_run
#     PROPERTIES
#     RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
# )




