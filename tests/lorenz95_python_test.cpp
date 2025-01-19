#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <gtest/gtest.h>
#include <Python.h>
#include <numpy/arrayobject.h>
#include <vector>
#include <stdexcept>

/**
 * @brief Initialize NumPy array API
 * @return true if initialization successful, false otherwise
 */
static bool init_numpy() {
    import_array();
    return true;
}

/**
 * @brief Test fixture for Lorenz95 Python integration tests
 * 
 * This fixture handles Python interpreter initialization/finalization and
 * provides access to the Lorenz95 Python module.
 */
class Lorenz95Test : public ::testing::Test {
protected:
    /** @brief Initialize Python interpreter and import Lorenz95 module */
    void SetUp() override {
        Py_Initialize();
        if (!init_numpy()) {
            throw std::runtime_error("NumPy initialization failed");
        }
        // Add the current directory to Python path
        PyRun_SimpleString("import sys; sys.path.append('.')");
        
        // Import the Lorenz95 module
        pModule = PyImport_ImportModule("lorenz95");
        if (!pModule) {
            PyErr_Print();
            throw std::runtime_error("Failed to import lorenz95 module");
        }
    }

    /** @brief Clean up Python resources */
    void TearDown() override {
        Py_XDECREF(pModule);
        Py_Finalize();
    }

    /** @brief Python module object for Lorenz95 */
    PyObject* pModule;
};

TEST_F(Lorenz95Test, TestLorenz95Integration) {
    // Get the integrate function
    PyObject* pFunc = PyObject_GetAttrString(pModule, "integrate");
    ASSERT_TRUE(pFunc != nullptr);

    // Create initial conditions (40 elements)
    std::vector<double> initial_state(40, 1.0);
    initial_state[19] = 1.001; // Small perturbation

    // Convert to NumPy array
    npy_intp dims[1] = {40};
    PyObject* pInitialState = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, initial_state.data());
    
    // Create arguments tuple
    PyObject* pArgs = PyTuple_New(3);
    PyTuple_SetItem(pArgs, 0, pInitialState);
    PyTuple_SetItem(pArgs, 1, PyFloat_FromDouble(0.05));  // dt
    PyTuple_SetItem(pArgs, 2, PyLong_FromLong(100));      // num_steps

    // Call the function
    PyObject* pResult = PyObject_CallObject(pFunc, pArgs);
    ASSERT_TRUE(pResult != nullptr);

    // Convert result back to C++
    PyArrayObject* np_arr = reinterpret_cast<PyArrayObject*>(pResult);
    double* result_data = static_cast<double*>(PyArray_DATA(np_arr));
    
    // Basic checks
    EXPECT_NE(result_data[0], initial_state[0]); // State should have changed
    EXPECT_FALSE(std::isnan(result_data[0]));    // Result should not be NaN

    // Cleanup
    Py_DECREF(pFunc);
    Py_DECREF(pArgs);
    Py_DECREF(pResult);
} 