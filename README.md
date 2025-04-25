# PreciseDotProduct

## Overview

PreciseDotProduct is a C++ project that focuses on computing the dot product of two vectors using precise floating point arithmetic. It leverages Boost.Multiprecision for high precision calculations and provides utilities for generating test cases, running tests, and verifying the accuracy of the dot product computation.

## Features

- Precise dot product calculation using custom algorithms.
- High precision dot product computation using Boost.Multiprecision.
- Test case generation for various vector characteristics.
- Test runner with detailed output and verification of results.

## Requirements

- Boost library

## Building the Project

1. Clone the repository:
   ```sh
   git clone https://github.com/eggunof/PreciseDotProduct.git
   cd PreciseDotProduct
   ```

2. Create a build directory and navigate to it:
   ```sh
   mkdir build && cd build
   ```

3. Run CMake to configure the project:
   ```sh
   cmake ..
   ```

4. Build the project:
   ```sh
   make
   ```

## Running the Tests

After building the project, the tests can be executed by running the generated executable:

```sh
./PreciseDotProduct
```

This will run the test suite and output the test results to the console.

## Code Structure

- `main.cc`: Contains the main logic for computing dot products, running tests, and utility functions.
- `CMakeLists.txt`: Configuration file for building the project with CMake.

## License

This project is licensed under the MIT License.