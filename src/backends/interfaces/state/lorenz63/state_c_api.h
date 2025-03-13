#pragma once

#ifdef __cplusplus
extern "C" {
#endif

// State functions
void* state_create(float x, float y, float z);
void state_destroy(void* state_ptr);
float state_get_x(void* state_ptr);
float state_get_y(void* state_ptr);
float state_get_z(void* state_ptr);
void state_set_x(void* state_ptr, float x);
void state_set_y(void* state_ptr, float y);
void state_set_z(void* state_ptr, float z);
float state_distance(void* state_ptr1, void* state_ptr2);

#ifdef __cplusplus
}
#endif