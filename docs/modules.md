# Using C++20 Modules in the Metada Framework

This document explains how to use the C++20 modules support in the framework, which offers several advantages over traditional header-based includes:

- Faster compilation times
- Better encapsulation and visibility control
- No macro leakage
- Clearer dependency management
- No need for header guards/include guards

## Modules Overview

The Metada Framework now offers C++20 module versions of the Geometry-related components:

- `metada.framework.interfaces.geometry` - The IGeometry interface
- `metada.framework.interfaces.geometry.iterator` - The IGeometryIterator interface
- `metada.framework.adapters.geometry` - The Geometry adapter implementation
- `metada.framework.adapters.geometry.iterator` - The GeometryIterator implementation
- `metada.framework.adapters.geometry.point_iterator` - The GeometryPointIterator implementation
- `metada.tests.mock.geometry` - The MockGeometry test implementation

## Using Modules in Your Code

To use C++20 modules in your own code, you can simply import the modules instead of including headers:

```cpp
// Old header approach:
#include "IGeometry.hpp"
#include "Geometry.hpp"

// New module approach:
import metada.framework.interfaces.geometry;
import metada.framework.adapters.geometry;
```

## Creating New Modules

To create a new module in the framework:

1. Create a .ixx file (or .cppm if you prefer) for your module interface
2. Use the `export module` syntax for the module declaration
3. Use `import` for dependencies
4. Use `export` to specify symbols visible outside the module
5. Update the relevant CMakeLists.txt to include your module file

Example module declaration:

```cpp
export module metada.your.module.name;

// Import other modules
import metada.framework.interfaces.geometry;

// Import standard library
import <vector>;
import <string>;

export namespace metada::your::namespace {
    // Exported symbols are visible to importers
    export class YourClass {
        // ...
    };
    
    // Non-exported symbols are only visible within the module
    class InternalHelper {
        // ...
    };
}
```

## CMake Configuration

To set up CMake for your module:

1. Add your module files to `target_sources`
2. Call `set_module_file_properties` with your module files
3. Call `add_module_compiler_flags` for the target

Example:

```cmake
target_sources(your_target
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/YourModule.ixx
)

set_module_file_properties("${CMAKE_CURRENT_SOURCE_DIR}/YourModule.ixx")
add_module_compiler_flags(your_target)
```

## Compiler Support

C++20 modules support varies by compiler:

- **GCC**: Limited support in GCC 10, better support in GCC 11+
- **Clang**: Limited support in Clang 13, better support in Clang 14+
- **MSVC**: Good support in VS 2019 16.10+ and VS 2022

The CMake configuration will warn you if your compiler's support for modules is limited.

## Transitioning Existing Code

When transitioning existing code to modules:

1. Create the corresponding .ixx file for your header
2. Move the implementation to the module file
3. Keep the original header for backward compatibility (with `#include` guards)
4. Gradually migrate client code to use `import` instead of `#include`

## Testing Module-based Code

The framework includes an updated `tests/framework/adapters/GeometryTest.cpp` that demonstrates testing module-based code.

## Known Limitations

- Debugging support may vary by IDE and compiler
- Build time improvements may not be significant for small projects
- Compatibility with some third-party libraries may be limited 