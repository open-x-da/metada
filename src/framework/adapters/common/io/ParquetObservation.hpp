#ifndef METADA_FRAMEWORK_PARQUET_OBSERVATION_HPP
#define METADA_FRAMEWORK_PARQUET_OBSERVATION_HPP

#include <memory>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

#include "DateTime.hpp"

// Forward declarations for Arrow
namespace arrow {
class Table;
class Column;
}  // namespace arrow

namespace metada::framework {

// Forward declaration for compatibility methods
struct ObsRecord;

/**
 * @brief Unified column-oriented observation data structure
 *
 * Replaces ObsRecord and handles Surface, Profile, and Satellite observations
 * with flexible schema support for different observation types.
 */
class ParquetObservation {
 public:
  // Observation types
  enum class ObsType {
    SURFACE,    // METAR, SYNOP, SHIPS, etc.
    PROFILE,    // Soundings, AIREP, etc.
    SATELLITE,  // Satellite radiances, retrievals
    UNKNOWN
  };

  // Data types for different fields
  using FieldValue = std::variant<double, int, std::string>;

  // Observation metadata
  struct ObservationMeta {
    std::string station_id;
    double longitude;
    double latitude;
    double elevation;
    double pressure;    // For profile/satellite, 0 for surface
    DateTime datetime;  // Using existing DateTime class
    int obs_type;       // WMO observation type code
    int channel;        // Satellite channel (0 for non-satellite)
    double obs_error;
    int qc_flag;
    std::string report_type;
    std::string instrument_type;
    ObsType category;  // Surface/Profile/Satellite
  };

  // Field definitions for different observation types
  struct FieldDefinition {
    std::string name;
    std::string unit;
    std::string description;
    bool required;
  };

  ParquetObservation();
  ~ParquetObservation() = default;

  // Core data access
  size_t size() const { return metadata_.size(); }
  bool empty() const { return metadata_.empty(); }

  // Metadata access
  const std::vector<ObservationMeta>& getMetadata() const { return metadata_; }
  const ObservationMeta& getMeta(size_t index) const;

  // Field access by name (returns vector of values for the field)
  std::vector<FieldValue> getField(const std::string& field_name) const;

  // Row access (for compatibility during transition)
  std::unordered_map<std::string, FieldValue> getObservation(
      size_t index) const;

  // Data manipulation
  void addObservation(
      const ObservationMeta& meta,
      const std::unordered_map<std::string, FieldValue>& fields);

  void clear();
  void reserve(size_t capacity);

  // Batch operations
  void appendObservations(const ParquetObservation& other);
  void filterByType(int obs_type);
  void filterByCategory(ObsType category);
  void filterByChannel(int channel);
  void filterByQC(int qc_flag);
  void filterByPressureRange(double min_pressure, double max_pressure);

  // Schema management
  void defineField(const std::string& field_name, const FieldDefinition& def);
  const std::unordered_map<std::string, FieldDefinition>& getSchema() const;

  // Parquet file operations
  bool saveToParquet(const std::string& filename) const;
  bool loadFromParquet(const std::string& filename);

  // Statistics and analysis
  std::vector<double> getPressureLevels() const;
  std::vector<int> getObservationTypes() const;
  size_t countByCategory(ObsType category) const;

  // Compatibility methods for transition
  std::vector<ObsRecord> toObsRecords() const;
  static ParquetObservation fromObsRecords(
      const std::vector<ObsRecord>& records);

 private:
  std::vector<ObservationMeta> metadata_;
  std::unordered_map<std::string, std::vector<FieldValue>> field_data_;
  std::unordered_map<std::string, FieldDefinition> schema_;

  // Helper methods
  void initializeDefaultSchema();
  void validateFieldData();
  void optimizeMemoryLayout();
  std::string getParquetSchema() const;
  void applyFilter(const std::vector<size_t>& indices_to_keep);
  void parseArrowTable(const std::shared_ptr<arrow::Table>& table);

  // Arrow/Parquet helper methods
  void extractDoubleColumn(const std::shared_ptr<arrow::Column>& column,
                           std::vector<double>& values);
  void extractIntColumn(const std::shared_ptr<arrow::Column>& column,
                        std::vector<int>& values);
  void extractStringColumn(const std::shared_ptr<arrow::Column>& column,
                           std::vector<std::string>& values);
  void extractFieldColumn(
      const std::shared_ptr<arrow::Column>& column,
      const std::string& field_name,
      std::unordered_map<std::string, std::vector<FieldValue>>& field_values);
};

}  // namespace metada::framework

#endif  // METADA_FRAMEWORK_PARQUET_OBSERVATION_HPP
