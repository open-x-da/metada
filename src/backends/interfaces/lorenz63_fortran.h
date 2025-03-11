#ifndef LORENZ63_FORTRAN_H
#define LORENZ63_FORTRAN_H

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

// Model functions
void* lorenz63_model_create(float sigma, float rho, float beta, float dt);
void lorenz63_model_destroy(void* model_ptr);
void lorenz63_model_get_params(void* model_ptr, float* sigma, float* rho, float* beta, float* dt);
void lorenz63_model_set_params(void* model_ptr, float sigma, float rho, float beta, float dt);

// Integrator functions
void* rk4_integrator_create(void* model_ptr);
void rk4_integrator_destroy(void* integrator_ptr);
void rk4_integrator_step(void* integrator_ptr, void* state_ptr);

#ifdef __cplusplus
}
#endif

#endif // LORENZ63_FORTRAN_H 