#include <string>
#include <chrono>
#include <stdexcept>

#include "ApplicationContext.hpp"
#include "AppTraits.hpp"
#include "../../../interfaces/Lorenz63.hpp"
#include "GoogleLogger.hpp"
#include "JsonConfig.hpp"

using namespace metada::backends::interfaces;
using namespace metada::framework::runs;
using namespace metada::backends::logger;
using namespace metada::backends::config;

using Traits = AppTraits<GoogleLogger, JsonConfig>;

int main() {
            // Create the application context
        auto app = ApplicationContext<Traits>("Lorenz63");

        // Get the logger and config
        auto& logger = app.getLogger();
        auto& config = app.getConfig();
    try {
        // Create a Lorenz63 state with initial values (1.0, 1.0, 1.0)
        Lorenz63State state(1.0f, 1.0f, 1.0f);
        Lorenz63State initial_state(1.0f, 1.0f, 1.0f);
        
        // Create the Lorenz63 model with standard parameters
        // sigma = 10.0, rho = 28.0, beta = 8/3, dt = 0.01
        Lorenz63Model model(10.0f, 28.0f, 8.0f/3.0f, 0.01f);
        
        // Get the model parameters to verify
        float sigma, rho, beta, dt;
        model.getParameters(sigma, rho, beta, dt);
        logger.Info() << "Model parameters:";
        logger.Info() << "  sigma = " << sigma;
        logger.Info() << "  rho = " << rho;
        logger.Info() << "  beta = " << beta;
        logger.Info() << "  dt = " << dt;
        
        // Create an integrator with the model
        Lorenz63Integrator integrator(model);
        
        // Measure performance
        auto start_time = std::chrono::high_resolution_clock::now();
        
        // Run the integration for 1000 steps
        const int num_steps = 1000;
        for (int i = 0; i < num_steps; ++i) {
            integrator.step(state);
            
            // Print every 100 steps
            if (i % 100 == 0) {
                logger.Info() << "Step " << i << ": " << state;
            }
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        // Calculate distance from initial to final state
        float distance = state.distance(initial_state);
        logger.Info() << "Distance from initial to final state: " << distance;
        
        // Print performance metrics
        logger.Info() << "Performance:";
        logger.Info() << "  Executed " << num_steps << " steps in " << duration.count() << " ms";
        logger.Info() << "  Average time per step: " << static_cast<float>(duration.count()) / num_steps << " ms";
        
        // Demonstrate vector access to components
        auto components = state.getComponents();
        logger.Info() << "Final state as vector: [" << components[0] << ", " << components[1] << ", " << components[2] << "]";
        
        // Final state using stream operator
        logger.Info() << "Final state using stream operator: " << state;
        
        // Demonstrate multiple states in one output
        logger.Info() << "Initial and final states: " << initial_state << " -> " << state;
    }
    catch (const std::exception& e) {
        logger.Error() << "Error: " << e.what();
        return 1;
    }
    
    return 0;
} 