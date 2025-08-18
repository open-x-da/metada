#include "ParquetObservation.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>

// Arrow/Parquet includes
#include <arrow/api.h>
#include <arrow/io/api.h>
#include <arrow/result.h>
#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>
#include <parquet/exception.h>

#include "ObsRecord.hpp"

namespace metada::framework {

ParquetObservation::ParquetObservation() {
  initializeDefaultSchema();
}

void ParquetObservation::initializeDefaultSchema() {
  // Surface observations
  defineField("temperature", {"temperature", "K", "Air temperature", true});
  defineField("humidity", {"humidity", "%", "Relative humidity", false});
  defineField("pressure", {"pressure", "hPa", "Surface pressure", true});
  defineField("wind_speed", {"wind_speed", "m/s", "Wind speed", false});
  defineField("wind_direction",
              {"wind_direction", "degrees", "Wind direction", false});

  // Profile observations
  defineField("height", {"height", "m", "Height above ground", false});
  defineField("geopotential",
              {"geopotential", "m²/s²", "Geopotential height", false});

  // Satellite observations
  defineField("brightness_temperature",
              {"brightness_temperature", "K", "Brightness temperature", false});
  defineField("radiance", {"radiance", "mW/m²/sr/cm⁻¹", "Radiance", false});
  defineField("cloud_fraction",
              {"cloud_fraction", "1", "Cloud fraction", false});
}

const ParquetObservation::ObservationMeta& ParquetObservation::getMeta(
    size_t index) const {
  if (index >= metadata_.size()) {
    throw std::out_of_range("Index out of range for metadata access");
  }
  return metadata_[index];
}

std::vector<ParquetObservation::FieldValue> ParquetObservation::getField(
    const std::string& field_name) const {
  auto it = field_data_.find(field_name);
  if (it == field_data_.end()) {
    return std::vector<FieldValue>();
  }
  return it->second;
}

std::unordered_map<std::string, ParquetObservation::FieldValue>
ParquetObservation::getObservation(size_t index) const {
  if (index >= metadata_.size()) {
    throw std::out_of_range("Index out of range for observation access");
  }

  std::unordered_map<std::string, FieldValue> result;

  // Add metadata fields
  const auto& meta = metadata_[index];
  result["station_id"] = meta.station_id;
  result["longitude"] = meta.longitude;
  result["latitude"] = meta.latitude;
  result["elevation"] = meta.elevation;
  result["pressure"] = meta.pressure;
  result["obs_type"] = meta.obs_type;
  result["channel"] = meta.channel;
  result["obs_error"] = meta.obs_error;
  result["qc_flag"] = meta.qc_flag;
  result["report_type"] = meta.report_type;
  result["instrument_type"] = meta.instrument_type;

  // Add field values
  for (const auto& [field_name, values] : field_data_) {
    if (index < values.size()) {
      result[field_name] = values[index];
    }
  }

  return result;
}

void ParquetObservation::addObservation(
    const ObservationMeta& meta,
    const std::unordered_map<std::string, FieldValue>& fields) {
  metadata_.push_back(meta);

  // Add field values, expanding vectors as needed
  for (const auto& [field_name, value] : fields) {
    if (field_data_.find(field_name) == field_data_.end()) {
      field_data_[field_name].resize(metadata_.size() - 1, FieldValue{});
    }
    field_data_[field_name].push_back(value);
  }

  // Fill missing fields with default values
  for (const auto& [field_name, _] : schema_) {
    if (field_data_[field_name].size() < metadata_.size()) {
      field_data_[field_name].resize(metadata_.size(), FieldValue{});
    }
  }
}

void ParquetObservation::clear() {
  metadata_.clear();
  field_data_.clear();
}

void ParquetObservation::reserve(size_t capacity) {
  metadata_.reserve(capacity);
  for (auto& [field_name, values] : field_data_) {
    values.reserve(capacity);
  }
}

void ParquetObservation::appendObservations(const ParquetObservation& other) {
  size_t old_size = metadata_.size();
  metadata_.insert(metadata_.end(), other.metadata_.begin(),
                   other.metadata_.end());

  // Append field data
  for (const auto& [field_name, values] : other.field_data_) {
    if (field_data_.find(field_name) == field_data_.end()) {
      field_data_[field_name].resize(old_size, FieldValue{});
    }
    field_data_[field_name].insert(field_data_[field_name].end(),
                                   values.begin(), values.end());
  }

  // Fill missing fields for existing observations
  for (auto& [field_name, values] : field_data_) {
    if (values.size() < metadata_.size()) {
      values.resize(metadata_.size(), FieldValue{});
    }
  }
}

void ParquetObservation::filterByType(int obs_type) {
  std::vector<size_t> indices_to_keep;
  for (size_t i = 0; i < metadata_.size(); ++i) {
    if (metadata_[i].obs_type == obs_type) {
      indices_to_keep.push_back(i);
    }
  }
  applyFilter(indices_to_keep);
}

void ParquetObservation::filterByCategory(ObsType category) {
  std::vector<size_t> indices_to_keep;
  for (size_t i = 0; i < metadata_.size(); ++i) {
    if (metadata_[i].category == category) {
      indices_to_keep.push_back(i);
    }
  }
  applyFilter(indices_to_keep);
}

void ParquetObservation::filterByChannel(int channel) {
  std::vector<size_t> indices_to_keep;
  for (size_t i = 0; i < metadata_.size(); ++i) {
    if (metadata_[i].channel == channel) {
      indices_to_keep.push_back(i);
    }
  }
  applyFilter(indices_to_keep);
}

void ParquetObservation::filterByQC(int qc_flag) {
  std::vector<size_t> indices_to_keep;
  for (size_t i = 0; i < metadata_.size(); ++i) {
    if (metadata_[i].qc_flag == qc_flag) {
      indices_to_keep.push_back(i);
    }
  }
  applyFilter(indices_to_keep);
}

void ParquetObservation::filterByPressureRange(double min_pressure,
                                               double max_pressure) {
  std::vector<size_t> indices_to_keep;
  for (size_t i = 0; i < metadata_.size(); ++i) {
    if (metadata_[i].pressure >= min_pressure &&
        metadata_[i].pressure <= max_pressure) {
      indices_to_keep.push_back(i);
    }
  }
  applyFilter(indices_to_keep);
}

void ParquetObservation::defineField(const std::string& field_name,
                                     const FieldDefinition& def) {
  schema_[field_name] = def;
}

const std::unordered_map<std::string, ParquetObservation::FieldDefinition>&
ParquetObservation::getSchema() const {
  return schema_;
}

bool ParquetObservation::saveToParquet(const std::string& filename) const {
  try {
    // Create Arrow arrays for each field
    std::vector<std::shared_ptr<arrow::Array>> arrays;
    std::vector<std::shared_ptr<arrow::Field>> fields;

    // Add metadata fields
    std::vector<double> lats, lons, pressures, times;
    std::vector<int> types, channels, qc_flags;
    std::vector<std::string> station_ids, report_types, instrument_types;
    std::vector<double> elevations, obs_errors;

    for (const auto& meta : metadata_) {
      lats.push_back(meta.latitude);
      lons.push_back(meta.longitude);
      pressures.push_back(meta.pressure);
      times.push_back(
          meta.datetime
              .toUnixTimestamp());  // Convert DateTime to Unix timestamp
      types.push_back(meta.obs_type);
      channels.push_back(meta.channel);
      qc_flags.push_back(meta.qc_flag);
      station_ids.push_back(meta.station_id);
      report_types.push_back(meta.report_type);
      instrument_types.push_back(meta.instrument_type);
      elevations.push_back(meta.elevation);
      obs_errors.push_back(meta.obs_error);
    }

    // Create Arrow arrays for metadata
    auto lat_array = std::make_shared<arrow::DoubleArray>(
        metadata_.size(),
        std::make_shared<arrow::DoubleBuffer>(lats.data(), lats.size()));
    auto lon_array = std::make_shared<arrow::DoubleArray>(
        metadata_.size(),
        std::make_shared<arrow::DoubleBuffer>(lons.data(), lons.size()));
    auto pressure_array = std::make_shared<arrow::DoubleArray>(
        metadata_.size(), std::make_shared<arrow::DoubleBuffer>(
                              pressures.data(), pressures.size()));
    auto time_array = std::make_shared<arrow::DoubleArray>(
        metadata_.size(),
        std::make_shared<arrow::DoubleBuffer>(times.data(), times.size()));
    auto type_array = std::make_shared<arrow::Int32Array>(
        metadata_.size(),
        std::make_shared<arrow::Int32Buffer>(types.data(), types.size()));
    auto channel_array = std::make_shared<arrow::Int32Array>(
        metadata_.size(),
        std::make_shared<arrow::Int32Buffer>(channels.data(), channels.size()));
    auto qc_array = std::make_shared<arrow::Int32Array>(
        metadata_.size(),
        std::make_shared<arrow::Int32Buffer>(qc_flags.data(), qc_flags.size()));
    auto elevation_array = std::make_shared<arrow::DoubleArray>(
        metadata_.size(), std::make_shared<arrow::DoubleBuffer>(
                              elevations.data(), elevations.size()));
    auto error_array = std::make_shared<arrow::DoubleArray>(
        metadata_.size(), std::make_shared<arrow::DoubleBuffer>(
                              obs_errors.data(), obs_errors.size()));

    // Create string arrays for metadata
    auto station_builder = std::make_shared<arrow::StringBuilder>();
    auto report_builder = std::make_shared<arrow::StringBuilder>();
    auto instrument_builder = std::make_shared<arrow::StringBuilder>();

    for (const auto& meta : metadata_) {
      PARQUET_THROW_NOT_OK(station_builder->Append(meta.station_id));
      PARQUET_THROW_NOT_OK(report_builder->Append(meta.report_type));
      PARQUET_THROW_NOT_OK(instrument_builder->Append(meta.instrument_type));
    }

    auto station_array = station_builder->Finish().ValueOrDie();
    auto report_array = report_builder->Finish().ValueOrDie();
    auto instrument_array = instrument_builder->Finish().ValueOrDie();

    // Add metadata fields to schema
    fields.push_back(
        std::make_shared<arrow::Field>("latitude", arrow::float64()));
    fields.push_back(
        std::make_shared<arrow::Field>("longitude", arrow::float64()));
    fields.push_back(
        std::make_shared<arrow::Field>("pressure", arrow::float64()));
    fields.push_back(std::make_shared<arrow::Field>("time", arrow::float64()));
    fields.push_back(
        std::make_shared<arrow::Field>("obs_type", arrow::int32()));
    fields.push_back(std::make_shared<arrow::Field>("channel", arrow::int32()));
    fields.push_back(std::make_shared<arrow::Field>("qc_flag", arrow::int32()));
    fields.push_back(
        std::make_shared<arrow::Field>("station_id", arrow::utf8()));
    fields.push_back(
        std::make_shared<arrow::Field>("report_type", arrow::utf8()));
    fields.push_back(
        std::make_shared<arrow::Field>("instrument_type", arrow::utf8()));
    fields.push_back(
        std::make_shared<arrow::Field>("elevation", arrow::float64()));
    fields.push_back(
        std::make_shared<arrow::Field>("obs_error", arrow::float64()));

    // Add metadata arrays
    arrays.push_back(lat_array);
    arrays.push_back(lon_array);
    arrays.push_back(pressure_array);
    arrays.push_back(time_array);
    arrays.push_back(type_array);
    arrays.push_back(channel_array);
    arrays.push_back(qc_array);
    arrays.push_back(station_array);
    arrays.push_back(report_array);
    arrays.push_back(instrument_array);
    arrays.push_back(elevation_array);
    arrays.push_back(error_array);

    // Create Arrow arrays for observation fields
    for (const auto& [field_name, values] : field_data_) {
      if (values.empty()) continue;

      // Handle different field types
      if (std::holds_alternative<double>(values[0])) {
        std::vector<double> double_values;
        for (const auto& val : values) {
          double_values.push_back(std::get<double>(val));
        }
        auto array = std::make_shared<arrow::DoubleArray>(
            values.size(), std::make_shared<arrow::DoubleBuffer>(
                               double_values.data(), double_values.size()));
        arrays.push_back(array);
        fields.push_back(
            std::make_shared<arrow::Field>(field_name, arrow::float64()));
      } else if (std::holds_alternative<int>(values[0])) {
        std::vector<int> int_values;
        for (const auto& val : values) {
          int_values.push_back(std::get<int>(val));
        }
        auto array = std::make_shared<arrow::Int32Array>(
            values.size(), std::make_shared<arrow::Int32Buffer>(
                               int_values.data(), int_values.size()));
        arrays.push_back(array);
        fields.push_back(
            std::make_shared<arrow::Field>(field_name, arrow::int32()));
      } else if (std::holds_alternative<std::string>(values[0])) {
        std::vector<std::string> string_values;
        for (const auto& val : values) {
          string_values.push_back(std::get<std::string>(val));
        }

        auto string_builder = std::make_shared<arrow::StringBuilder>();
        for (const auto& str_val : string_values) {
          PARQUET_THROW_NOT_OK(string_builder->Append(str_val));
        }

        auto array = string_builder->Finish().ValueOrDie();
        arrays.push_back(array);
        fields.push_back(
            std::make_shared<arrow::Field>(field_name, arrow::utf8()));
      }
    }

    // Create Arrow table
    auto schema = std::make_shared<arrow::Schema>(fields);
    auto table = arrow::Table::Make(schema, arrays);

    // Write to Parquet file
    std::shared_ptr<arrow::io::FileOutputStream> outfile;
    PARQUET_ASSIGN_OR_THROW(outfile,
                            arrow::io::FileOutputStream::Open(filename));

    PARQUET_THROW_NOT_OK(parquet::arrow::WriteTable(
        *table, arrow::default_memory_pool(), outfile));

    std::cout << "Successfully saved " << size()
              << " observations to Parquet file: " << filename << std::endl;
    return true;
  } catch (const std::exception& e) {
    std::cerr << "Error saving to Parquet: " << e.what() << std::endl;
    return false;
  }
}

bool ParquetObservation::loadFromParquet(const std::string& filename) {
  try {
    // Read from Parquet file using Arrow
    std::shared_ptr<arrow::io::ReadableFile> infile;
    PARQUET_ASSIGN_OR_THROW(infile, arrow::io::ReadableFile::Open(filename));

    std::unique_ptr<parquet::arrow::FileReader> reader;
    PARQUET_THROW_NOT_OK(parquet::arrow::OpenFile(
        infile, arrow::default_memory_pool(), &reader));

    std::shared_ptr<arrow::Table> table;
    PARQUET_THROW_NOT_OK(reader->ReadTable(&table));

    // Parse table back to internal structure
    parseArrowTable(table);

    std::cout << "Successfully loaded " << size()
              << " observations from Parquet file: " << filename << std::endl;
    return true;
  } catch (const std::exception& e) {
    std::cerr << "Error loading from Parquet: " << e.what() << std::endl;
    return false;
  }
}

void ParquetObservation::parseArrowTable(
    const std::shared_ptr<arrow::Table>& table) {
  // Clear existing data
  clear();

  // Get table dimensions
  int64_t num_rows = table->num_rows();
  int num_cols = table->num_columns();

  if (num_rows == 0) return;

  // Reserve space for efficiency
  reserve(num_rows);

  // Extract metadata columns
  std::vector<double> lats, lons, pressures, times;
  std::vector<int> types, channels, qc_flags;
  std::vector<std::string> station_ids, report_types, instrument_types;
  std::vector<double> elevations, obs_errors;

  // Extract observation field data
  std::unordered_map<std::string, std::vector<FieldValue>> field_values;

  // Process each column
  for (int col = 0; col < num_cols; ++col) {
    auto column = table->column(col);
    auto field_name = table->schema()->field(col)->name();

    if (field_name == "latitude") {
      extractDoubleColumn(column, lats);
    } else if (field_name == "longitude") {
      extractDoubleColumn(column, lons);
    } else if (field_name == "pressure") {
      extractDoubleColumn(column, pressures);
    } else if (field_name == "time") {
      extractDoubleColumn(column, times);
    } else if (field_name == "obs_type") {
      extractIntColumn(column, types);
    } else if (field_name == "channel") {
      extractIntColumn(column, channels);
    } else if (field_name == "qc_flag") {
      extractIntColumn(column, qc_flags);
    } else if (field_name == "station_id") {
      extractStringColumn(column, station_ids);
    } else if (field_name == "report_type") {
      extractStringColumn(column, report_types);
    } else if (field_name == "instrument_type") {
      extractStringColumn(column, instrument_types);
    } else if (field_name == "elevation") {
      extractDoubleColumn(column, elevations);
    } else if (field_name == "obs_error") {
      extractDoubleColumn(column, obs_errors);
    } else {
      // This is an observation field
      extractFieldColumn(column, field_name, field_values);
    }
  }

  // Create ObservationMeta objects and populate field_data_
  for (int64_t row = 0; row < num_rows; ++row) {
    ObservationMeta meta;
    meta.latitude = lats[row];
    meta.longitude = lons[row];
    meta.pressure = pressures[row];
    meta.datetime = DateTime::fromUnixTimestamp(times[row]);
    meta.obs_type = types[row];
    meta.channel = channels[row];
    meta.qc_flag = qc_flags[row];
    meta.station_id = station_ids[row];
    meta.report_type = report_types[row];
    meta.instrument_type = instrument_types[row];
    meta.elevation = elevations[row];
    meta.obs_error = obs_errors[row];

    // Determine category based on pressure
    if (meta.pressure > 0) {
      meta.category = ObsType::PROFILE;
    } else {
      meta.category = ObsType::SURFACE;
    }

    metadata_.push_back(meta);

    // Add field values for this observation
    for (const auto& [field_name, values] : field_values) {
      if (row < values.size()) {
        field_data_[field_name].push_back(values[row]);
      }
    }
  }
}

void ParquetObservation::extractDoubleColumn(
    const std::shared_ptr<arrow::Column>& column, std::vector<double>& values) {
  auto array =
      std::static_pointer_cast<arrow::DoubleArray>(column->data()->chunk(0));
  values.clear();
  values.reserve(array->length());

  for (int64_t i = 0; i < array->length(); ++i) {
    if (array->IsValid(i)) {
      values.push_back(array->Value(i));
    } else {
      values.push_back(0.0);  // Default value for missing data
    }
  }
}

void ParquetObservation::extractIntColumn(
    const std::shared_ptr<arrow::Column>& column, std::vector<int>& values) {
  auto array =
      std::static_pointer_cast<arrow::Int32Array>(column->data()->chunk(0));
  values.clear();
  values.reserve(array->length());

  for (int64_t i = 0; i < array->length(); ++i) {
    if (array->IsValid(i)) {
      values.push_back(array->Value(i));
    } else {
      values.push_back(0);  // Default value for missing data
    }
  }
}

void ParquetObservation::extractStringColumn(
    const std::shared_ptr<arrow::Column>& column,
    std::vector<std::string>& values) {
  auto array =
      std::static_pointer_cast<arrow::StringArray>(column->data()->chunk(0));
  values.clear();
  values.reserve(array->length());

  for (int64_t i = 0; i < array->length(); ++i) {
    if (array->IsValid(i)) {
      values.push_back(array->GetString(i));
    } else {
      values.push_back("");  // Default value for missing data
    }
  }
}

void ParquetObservation::extractFieldColumn(
    const std::shared_ptr<arrow::Column>& column, const std::string& field_name,
    std::unordered_map<std::string, std::vector<FieldValue>>& field_values) {
  auto data_type = column->type();
  std::vector<FieldValue> values;

  if (data_type->id() == arrow::Type::DOUBLE) {
    auto array =
        std::static_pointer_cast<arrow::DoubleArray>(column->data()->chunk(0));
    values.reserve(array->length());

    for (int64_t i = 0; i < array->length(); ++i) {
      if (array->IsValid(i)) {
        values.push_back(array->Value(i));
      } else {
        values.push_back(0.0);
      }
    }
  } else if (data_type->id() == arrow::Type::INT32) {
    auto array =
        std::static_pointer_cast<arrow::Int32Array>(column->data()->chunk(0));
    values.reserve(array->length());

    for (int64_t i = 0; i < array->length(); ++i) {
      if (array->IsValid(i)) {
        values.push_back(array->Value(i));
      } else {
        values.push_back(0);
      }
    }
  } else if (data_type->id() == arrow::Type::STRING) {
    auto array =
        std::static_pointer_cast<arrow::StringArray>(column->data()->chunk(0));
    values.reserve(array->length());

    for (int64_t i = 0; i < array->length(); ++i) {
      if (array->IsValid(i)) {
        values.push_back(array->GetString(i));
      } else {
        values.push_back("");
      }
    }
  }

  if (!values.empty()) {
    field_values[field_name] = std::move(values);
  }
}

std::vector<double> ParquetObservation::getPressureLevels() const {
  std::vector<double> levels;
  levels.reserve(metadata_.size());

  for (const auto& meta : metadata_) {
    if (meta.pressure > 0) {  // Only include valid pressure levels
      levels.push_back(meta.pressure);
    }
  }

  // Remove duplicates and sort
  std::sort(levels.begin(), levels.end());
  levels.erase(std::unique(levels.begin(), levels.end()), levels.end());

  return levels;
}

std::vector<int> ParquetObservation::getObservationTypes() const {
  std::vector<int> types;
  types.reserve(metadata_.size());

  for (const auto& meta : metadata_) {
    types.push_back(meta.obs_type);
  }

  // Remove duplicates and sort
  std::sort(types.begin(), types.end());
  types.erase(std::unique(types.begin(), types.end()), types.end());

  return types;
}

size_t ParquetObservation::countByCategory(ObsType category) const {
  size_t count = 0;
  for (const auto& meta : metadata_) {
    if (meta.category == category) {
      ++count;
    }
  }
  return count;
}

std::vector<ObsRecord> ParquetObservation::toObsRecords() const {
  // TODO: Implement conversion to ObsRecord format for compatibility
  // This will be needed during the transition period
  std::vector<ObsRecord> records;
  // Implementation will go here
  return records;
}

ParquetObservation ParquetObservation::fromObsRecords(
    const std::vector<ObsRecord>& records) {
  // TODO: Implement conversion from ObsRecord format
  // This will be needed during the transition period
  ParquetObservation result;
  // Implementation will go here
  return result;
}

void ParquetObservation::validateFieldData() {
  // TODO: Implement field data validation
  // Check that all required fields are present and have correct types
}

void ParquetObservation::optimizeMemoryLayout() {
  // TODO: Implement memory layout optimization
  // This could include reordering columns for better cache locality
}

std::string ParquetObservation::getParquetSchema() const {
  // TODO: Implement Parquet schema generation
  // This will be used when writing to Parquet files
  return "parquet_schema_placeholder";
}

// Private helper method for filtering
void ParquetObservation::applyFilter(
    const std::vector<size_t>& indices_to_keep) {
  if (indices_to_keep.empty()) {
    clear();
    return;
  }

  // Create new metadata vector
  std::vector<ObservationMeta> new_metadata;
  new_metadata.reserve(indices_to_keep.size());
  for (size_t idx : indices_to_keep) {
    new_metadata.push_back(metadata_[idx]);
  }
  metadata_ = std::move(new_metadata);

  // Update field data
  for (auto& [field_name, values] : field_data_) {
    std::vector<FieldValue> new_values;
    new_values.reserve(indices_to_keep.size());
    for (size_t idx : indices_to_keep) {
      if (idx < values.size()) {
        new_values.push_back(values[idx]);
      } else {
        new_values.push_back(FieldValue{});
      }
    }
    values = std::move(new_values);
  }
}

}  // namespace metada::framework
