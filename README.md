## Building the Project

### Dependencies

- **Intel oneAPI Base Toolkit & Intel oneAPI HPC Toolkit**: Provides Intel IPP (Integrated Performance Primitives) and MKL (Math Kernel Library) for FFT operations and optimized vector functions.
- **CMake** (version 3.12 or higher): For cross-platform build configuration.
- **OpenMP**: For shared-memory parallelization (usually included with modern compilers).
- **C++17 compliant compiler**: e.g., GCC 7+, Clang 5+, MSVC 2017+.

### Build Instructions

#### Using CMake (Recommended for all platforms)

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/yourusername/pdb2mrc.git
    cd pdb2mrc
    ```

2.  **Create a build directory:**
    ```bash
    mkdir build
    cd build
    ```

3.  **Configure the project with CMake:**
    - On **Linux/macOS**:
        Ensure the Intel oneAPI environment is sourced. The CMake script is designed to find Intel IPP and MKL automatically if they are in standard paths.
        ```bash
        source /opt/intel/oneapi/setvars.sh
        cmake .. -DCMAKE_BUILD_TYPE=Release
        ```
    - On **Windows**:
        Open "Intel oneAPI command prompt for x64" and run:
        ```bash
        cmake .. -G "Visual Studio 17 2022" -A x64
        ```
        Or for Ninja:
        ```bash
        cmake .. -G Ninja -DCMAKE_BUILD_TYPE=Release
        ```

4.  **Build the project:**
    ```bash
    cmake --build . --config Release
    ```

5.  **Install (optional):**
    ```bash
    cmake --install . --prefix /path/to/install
    ```

#### Using Visual Studio (Windows only)

1.  Open the "Intel oneAPI command prompt for x64".
2.  Navigate to the project root.
3.  Generate Visual Studio solution files using CMake as described above.
4.  Open the generated `pdb2mrc.sln` file in Visual Studio and build.

#### Using Makefiles (Linux/macOS)

After configuring with CMake, you can simply run `make` in the build directory.
