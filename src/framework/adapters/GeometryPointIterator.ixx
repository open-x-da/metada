/**
 * @file GeometryPointIterator.ixx
 * @brief Point iterator implementation for Geometry adapter (module version)
 * @ingroup adapters
 *
 * @details
 * This module extends the GeometryIterator to provide specific functionality
 * for iterating over grid points in a Geometry, with proper coordinate updates.
 */

export module metada.framework.adapters.geometry.point_iterator;

// Import dependencies
import metada.framework.adapters.geometry.iterator;
import metada.framework.interfaces.geometry.iterator;

// Standard library includes (using traditional #include instead of import)
#include <functional>
#include <vector>

// Forward declaration of Geometry class
export namespace metada::framework {

template <typename Backend>
class Geometry;

/**
 * @brief Extended iterator implementation for traversing geometry grid points
 * with coordinate updates
 *
 * @tparam T Type of the coordinate value (typically double)
 */
template <typename T>
class GeometryPointIterator : public GeometryIterator<T> {
 public:
  /**
   * @brief Default constructor creates an end iterator
   */
  GeometryPointIterator() = default;

  /**
   * @brief Construct iterator with specific position and geometry
   *
   * @param position Current position indices
   * @param dimensions Grid dimensions
   * @param coordinates Initial coordinates
   * @param geometry Pointer to the geometry object
   */
  template <typename Backend>
  GeometryPointIterator(const std::vector<size_t>& position,
                        const std::vector<size_t>& dimensions,
                        const std::vector<T>& coordinates,
                        const Geometry<Backend>* geometry)
      : GeometryIterator<T>(position, dimensions, coordinates),
        geometry_(const_cast<Geometry<Backend>*>(geometry)),
        get_coordinates_([this](const std::vector<size_t>& pos) {
          return static_cast<Geometry<Backend>*>(geometry_)->getCoordinates(
              pos);
        }) {}

  /**
   * @brief Copy constructor
   */
  GeometryPointIterator(const GeometryPointIterator& other)
      : GeometryIterator<T>(other),
        geometry_(other.geometry_),
        get_coordinates_(other.get_coordinates_) {}

  /**
   * @brief Pre-increment operator with coordinate update
   */
  IGeometryIterator<T>& operator++() override {
    // Using the base class operator++
    GeometryIterator<T>::operator++();

    // Update coordinates after incrementing position
    if (geometry_ && !GeometryIterator<T>::isDone()) {
      auto& pos = GeometryIterator<T>::getPosition();
      auto coords = get_coordinates_(pos);
      GeometryIterator<T>::coordinates_ = coords;
    }

    return *this;
  }

  /**
   * @brief Clone this iterator
   * @return Unique pointer to a new iterator instance
   */
  std::unique_ptr<IGeometryIterator<T>> clone() const override {
    return std::make_unique<GeometryPointIterator<T>>(*this);
  }

 private:
  void* geometry_ = nullptr;  ///< Generic pointer to any Geometry<Backend>
  std::function<std::vector<T>(const std::vector<size_t>&)> get_coordinates_;
};

}  // namespace metada::framework 