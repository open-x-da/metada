#pragma once

#include <iterator>
#include <vector>

#include "PointObservation.hpp"

namespace metada::backends::common::observation {

using ObservationPoint = metada::framework::ObservationPoint;

class PrepBUFRObservationIterator {
 public:
  using iterator_category = std::forward_iterator_tag;
  using value_type = ObservationPoint;
  using difference_type = std::ptrdiff_t;
  using pointer = ObservationPoint*;
  using reference = ObservationPoint&;

  PrepBUFRObservationIterator() = default;

  PrepBUFRObservationIterator(const std::vector<ObservationPoint>* data,
                              size_t index)
      : data_(data), index_(index) {}

  reference operator*() { return const_cast<reference>((*data_)[index_]); }
  pointer operator->() { return &const_cast<reference>((*data_)[index_]); }

  PrepBUFRObservationIterator& operator++() {
    ++index_;
    return *this;
  }
  PrepBUFRObservationIterator operator++(int) {
    PrepBUFRObservationIterator tmp = *this;
    ++index_;
    return tmp;
  }

  bool operator==(const PrepBUFRObservationIterator& other) const {
    return data_ == other.data_ && index_ == other.index_;
  }
  bool operator!=(const PrepBUFRObservationIterator& other) const {
    return !(*this == other);
  }

 private:
  const std::vector<ObservationPoint>* data_{nullptr};
  size_t index_{0};
};

}  // namespace metada::backends::common::observation
