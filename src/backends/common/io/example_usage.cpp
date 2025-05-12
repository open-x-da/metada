#include <iostream>
#include <string>

#include "BufrFortranWrapper.hpp"

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <bufr_filename>" << std::endl;
    return 1;
  }

  std::string filename = argv[1];

  try {
    // Method 1: Let the wrapper allocate a unit number
    {
      std::cout << "Method 1: Using automatic unit number allocation"
                << std::endl;
      metada::backends::io::BufrFortranWrapper bufr;
      bufr.open(filename);

      std::string subset;
      int date;

      std::cout << "Using unit number: " << bufr.getUnitNumber() << std::endl;

      // Read and process data
      while (bufr.readNextSubset(subset, date)) {
        std::cout << "Read subset: " << subset << ", date: " << date
                  << std::endl;
      }

      // File will be automatically closed when bufr goes out of scope
    }

    // Method 2: Specify a unit number explicitly (Fortran-style)
    {
      std::cout << "\nMethod 2: Using explicit unit number (Fortran-style)"
                << std::endl;
      metada::backends::io::BufrFortranWrapper bufr;

      // Open with specific unit number 11 (like in Fortran)
      bufr.openWithUnit(11, filename);

      std::string subset;
      int date;

      std::cout << "Using unit number: " << bufr.getUnitNumber() << std::endl;

      // Read and process data
      while (bufr.readNextSubset(subset, date)) {
        std::cout << "Read subset: " << subset << ", date: " << date
                  << std::endl;
      }

      // File will be automatically closed when bufr goes out of scope
    }

    std::cout << "Successfully processed file: " << filename << std::endl;
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
}