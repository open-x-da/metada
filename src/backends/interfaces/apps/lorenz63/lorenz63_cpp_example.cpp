#include <iostream>
#include <iomanip>
#include "../../../interfaces/Lorenz63.hpp"

using namespace metada::backends::interfaces;

int main() {
    // Create a Lorenz63 state with initial values (1.0, 1.0, 1.0)
    Lorenz63State state(1.0f, 1.0f, 1.0f);
    Lorenz63State initial_state(1.0f, 1.0f, 1.0f);
    
    // Print initial state
    std::cout << "Initial state:" << std::endl;
    std::cout << "  x = " << state.getX() << std::endl;
    std::cout << "  y = " << state.getY() << std::endl;
    std::cout << "  z = " << state.getZ() << std::endl;
    
    // Create the Lorenz63 model with standard parameters
    // sigma = 10.0, rho = 28.0, beta = 8/3, dt = 0.01
    Lorenz63Model model(10.0f, 28.0f, 8.0f/3.0f, 0.01f);
    
    // Get the model parameters to verify
    float sigma, rho, beta, dt;
    model.getParameters(sigma, rho, beta, dt);
    std::cout << "\nModel parameters:" << std::endl;
    std::cout << "  sigma = " << sigma << std::endl;
    std::cout << "  rho = " << rho << std::endl;
    std::cout << "  beta = " << beta << std::endl;
    std::cout << "  dt = " << dt << std::endl;
    
    // Create an integrator with the model
    Lorenz63Integrator integrator(model);
    
    // Run the integration for 1000 steps
    const int num_steps = 1000;
    for (int i = 0; i < num_steps; ++i) {
        integrator.step(state);
        
        // Print every 100 steps
        if (i % 100 == 0) {
            std::cout << "\nStep " << i << ":" << std::endl;
            std::cout << "  x = " << state.getX() << std::endl;
            std::cout << "  y = " << state.getY() << std::endl;
            std::cout << "  z = " << state.getZ() << std::endl;
        }
    }
    
    // Calculate distance from initial to final state
    float distance = state.distance(initial_state);
    std::cout << "\nDistance from initial to final state: " << distance << std::endl;
    
    return 0;
} 