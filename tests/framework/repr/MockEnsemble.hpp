#ifndef METADA_TESTS_FRAMEWORK_REPR_MOCK_ENSEMBLE_HPP_
#define METADA_TESTS_FRAMEWORK_REPR_MOCK_ENSEMBLE_HPP_

#include <gmock/gmock.h>

#include "IEnsemble.hpp"

namespace metada {
namespace framework {
namespace repr {
namespace tests {

class MockEnsemble : public IEnsemble {
 public:
  // Ensemble management
  MOCK_METHOD(void, initialize, (const tools::config::IConfig& config),
              (override));
  MOCK_METHOD(size_t, getSize, (), (const, override));
  MOCK_METHOD(void, resize, (size_t new_size), (override));

  // Member access
  MOCK_METHOD(IState&, getMember, (size_t index), (override));
  MOCK_METHOD(const IState&, getMember, (size_t index), (const, override));

  // Statistical operations
  MOCK_METHOD(void, computeMean, (), (override));
  MOCK_METHOD(const IState&, getMean, (), (const, override));
  MOCK_METHOD(void, computePerturbations, (), (override));
  MOCK_METHOD(const IState&, getPerturbation, (size_t index),
              (const, override));

  // LETKF specific operations
  MOCK_METHOD(void, inflate, (double factor), (override));
  MOCK_METHOD(void, transform,
              (const std::vector<std::vector<double>>& transform_matrix),
              (override));
  MOCK_METHOD(void, localizeCovariance,
              (const std::vector<double>& localization_weights), (override));

  // Validation
  MOCK_METHOD(void, validate, (), (const, override));
  MOCK_METHOD(bool, isValid, (), (const, override));
};

}  // namespace tests
}  // namespace repr
}  // namespace framework
}  // namespace metada

#endif  // METADA_TESTS_FRAMEWORK_REPR_MOCK_ENSEMBLE_HPP_