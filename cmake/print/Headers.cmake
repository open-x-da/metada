# Print section headers
function(print_section_header text)
    message(STATUS "")
    message(STATUS "${text}")
    string(LENGTH "${text}" text_length)
    string(REPEAT "-" ${text_length} separator)
    message(STATUS "${separator}")
    message(STATUS "")
endfunction() 