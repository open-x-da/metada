# Define the scientific representation library
add_library(repr INTERFACE
    State.hpp
    Model.hpp
    Observation.hpp
    Increment.hpp
    ObsOperator.hpp
)

# Set up include directories
target_include_directories(repr
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
)

# Create alias target
add_library(metada::framework::repr ALIAS repr)

# Install headers
install(
    FILES
        State.hpp
        Model.hpp
        Observation.hpp
        Increment.hpp
        ObsOperator.hpp
    DESTINATION
        include/metada/framework/repr
)

# Install library
install(
    TARGETS repr
    EXPORT metada-targets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
) 