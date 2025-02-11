# Define the scientific representation library
add_library(repr_interface INTERFACE)

# Set up include directories
target_include_directories(repr_interface
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
)

# Create alias target
add_library(metada::framework::repr ALIAS repr_interface)

# Install headers
install(
    FILES
        IState.hpp
        State.hpp
        IModel.hpp
        Model.hpp
        IObservation.hpp
        Observation.hpp
        IIncrement.hpp
        Increment.hpp
        IObsOperator.hpp
        ObsOperator.hpp
    DESTINATION
        include/metada/framework/repr
)

# Install library
install(
    TARGETS repr_interface
    EXPORT metada-targets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
) 