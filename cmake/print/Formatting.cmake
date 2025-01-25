# Formatting utilities for print functions
function(print_header text)
    message(STATUS "\n=== ${text} ===\n")
endfunction()

function(print_subheader text)
    message(STATUS "\n--- ${text} ---")
endfunction()

function(print_item text)
    message(STATUS "  ${text}")
endfunction()

function(print_subitem text)
    message(STATUS "    ${text}")
endfunction()