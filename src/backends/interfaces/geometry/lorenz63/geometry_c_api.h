#pragma once

#ifdef __cplusplus
extern "C" {
#endif

// Geometry functions
void* geometry_create(float x_min, float x_max, float y_min, float y_max, float z_min, float z_max);
void geometry_destroy(void* geometry_ptr);
int geometry_contains_point(void* geometry_ptr, void* state_ptr);
void geometry_get_x_range(void* geometry_ptr, float* min_val, float* max_val);
void geometry_get_y_range(void* geometry_ptr, float* min_val, float* max_val);
void geometry_get_z_range(void* geometry_ptr, float* min_val, float* max_val);

#ifdef __cplusplus
}
#endif 