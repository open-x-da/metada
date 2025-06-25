#include "SimpleStateIterator.hpp"

#include "SimpleState.hpp"

namespace metada::backends::simple {

SimpleStateIterator::reference SimpleStateIterator::operator*() {
  const auto& coord = *geom_it_;
  auto grid_coord = coord.getGridCoords2D();
  current_value_ = std::make_pair(coord, state_->at(grid_coord));
  return current_value_;
}

}  // namespace metada::backends::simple