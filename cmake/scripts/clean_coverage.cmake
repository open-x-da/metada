# Get all .gcda files recursively
file(GLOB_RECURSE coverage_files
  "${CMAKE_CURRENT_BINARY_DIR}/*.gcda"
)

# Remove each coverage file
foreach(file ${coverage_files})
  file(REMOVE ${file})
endforeach() 