/**
 * @file PointObservation.hpp
 * @brief Generic data structures for point-based observations
 * @ingroup adapters
 * @author Metada Framework Team
 *
 * @details
 * This file contains generic data structures for representing point-based
 * observations in meteorological and oceanographic data assimilation systems.
 * The structures include spatial location information, observation values,
 * associated errors, and validity flags. These components are designed to be
 * reused across different observation-related backend implementations to
 * promote consistency, interoperability, and code reuse throughout the
 * framework.
 *
 * @note These structures are designed to work seamlessly with the Location
 *       class and various coordinate systems supported by the framework.
 */
#pragma once

#include "Location.hpp"

namespace metada::framework {

/**
 * @brief Point-based observation data structure with location and measurement
 * information
 *
 * The ObservationPoint class encapsulates all information associated with a
 * single point observation, including its spatial location, measured value,
 * measurement uncertainty, and data quality flag. This class serves as a
 * fundamental building block for observation operators and data assimilation
 * algorithms.
 *
 * The class is designed to be lightweight and efficient, suitable for handling
 * large volumes of observational data in real-time and batch processing
 * scenarios. It integrates seamlessly with the Location class to support
 * multiple coordinate systems (grid, geographic, and Cartesian).
 *
 * @example
 * @code
 * // Create observation at geographic location
 * Location geoLoc(45.0, -122.0, 1000.0, CoordinateSystem::GEOGRAPHIC);
 * ObservationPoint tempObs(geoLoc, 15.5, 0.2);  // 15.5°C ± 0.2°C
 *
 * // Create observation at grid location
 * Location gridLoc(10, 20, 5);
 * ObservationPoint pressObs(gridLoc, 1013.25, 1.0);  // 1013.25 hPa ± 1.0 hPa
 *
 * // Create invalid observation placeholder
 * ObservationPoint invalidObs(geoLoc);  // Invalid observation at location
 * @endcode
 *
 * @note Member variables are kept public for direct access and compatibility
 *       with existing serialization and data processing workflows.
 */
class ObservationPoint {
 public:
  // ============================================================================
  // MEMBER VARIABLES
  // ============================================================================

  Location location;  ///< Spatial location of the observation using framework's
                      ///< Location class
  double value;       ///< Measured/observed value in appropriate physical units
  double error;       ///< Observation error or uncertainty (standard deviation)
  bool is_valid;  ///< Data quality flag indicating whether the observation is
                  ///< valid for use

  // ============================================================================
  // CONSTRUCTORS
  // ============================================================================

  /**
   * @brief Construct a valid observation with location, value, and error
   *
   * Creates a complete observation point with all necessary information for
   * data assimilation. The observation is automatically marked as valid.
   *
   * @param loc Spatial location of the observation
   * @param val Measured value in appropriate physical units
   * @param err Observation error/uncertainty (typically standard deviation)
   *
   * @note The error parameter should represent the standard deviation of the
   *       measurement uncertainty, not variance or other error metrics.
   */
  ObservationPoint(const Location& loc, double val, double err)
      : location(loc), value(val), error(err), is_valid(true) {}

  /**
   * @brief Construct an invalid observation placeholder at a location
   *
   * Creates an observation point with only location information, marking it
   * as invalid. This constructor is useful for pre-allocating observation
   * structures or representing missing/rejected observations.
   *
   * @param loc Spatial location where observation should be taken
   *
   * @note The value and error are initialized to 0.0, and is_valid is set to
   * false. These should be set appropriately before using the observation.
   */
  ObservationPoint(const Location& loc)
      : location(loc), value(0.0), error(0.0), is_valid(false) {}

  // ============================================================================
  // COMPARISON OPERATORS
  // ============================================================================

  /**
   * @brief Equality comparison operator for ObservationPoint objects
   *
   * Compares all components of two observation points for exact equality:
   * location, value, error, and validity flag. This is useful for testing,
   * validation, and duplicate detection in observation datasets.
   *
   * @param other ObservationPoint to compare with
   * @return true if all components (location, value, error, is_valid) are
   * identical
   * @return false if any component differs
   *
   * @note This performs exact floating-point comparison for value and error.
   *       Consider using approximate comparison for real-world applications
   *       where floating-point precision may be an issue.
   */
  bool operator==(const ObservationPoint& other) const {
    return location == other.location && value == other.value &&
           error == other.error && is_valid == other.is_valid;
  }
};

}  // namespace metada::framework