import numpy as np

def lorenz95(x, dt):
    ##
    # @brief Single step of Lorenz 95 model
    # @param x Input state array
    # @param dt Time step
    # @return Updated state array
    #
    N = len(x)
    d = np.zeros(N, dtype=np.float64)
    F = 8.0
    
    # Compute derivatives using vectorized operations
    # This helps prevent overflow by better handling the calculations
    x_plus1 = np.roll(x, -1)  # x[i+1]
    x_minus2 = np.roll(x, 2)  # x[i-2]
    x_minus1 = np.roll(x, 1)  # x[i-1]
    
    # Compute the tendency term with better numerical stability
    d = (x_plus1 - x_minus2) * x_minus1 - x + F
    
    # Apply limiting to prevent extreme values
    max_val = 1e10
    d = np.clip(d, -max_val, max_val)
    
    # Use more stable integration
    new_x = x + dt * d
    
    # Clip the result to prevent unbounded growth
    return np.clip(new_x, -max_val, max_val)

def integrate(initial_state, dt, num_steps):
    """Integrate Lorenz 95 model for multiple steps"""
    state = np.array(initial_state, dtype=np.float64)
    for _ in range(num_steps):
        state = lorenz95(state, dt)
    return state 