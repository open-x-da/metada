#pragma once

#ifdef __cplusplus
extern "C" {
#endif

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