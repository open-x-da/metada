# 递归查找所有 .gcda 文件
file(GLOB_RECURSE GCDA_FILES "${CMAKE_BINARY_DIR}/*.gcda")

# 输出找到的文件数量
message(STATUS "Found ${GCDA_FILES} .gcda files to delete.")

# 删除所有 .gcda 文件
if(GCDA_FILES)
    file(REMOVE ${GCDA_FILES})
    message(STATUS "Deleted all .gcda files.")
else()
    message(STATUS "No .gcda files found.")
endif()