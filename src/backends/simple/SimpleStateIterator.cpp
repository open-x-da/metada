#include "SimpleStateIterator.hpp"

#include "SimpleState.hpp"

namespace metada::backends::simple {

SimpleStateIterator::reference SimpleStateIterator::operator*() {
  const auto& coord = *geom_it_;
  current_value_ = std::make_pair(coord, state_->at(coord));
  return current_value_;
}

}  // namespace metada::backends::simple