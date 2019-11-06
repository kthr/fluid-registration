Dependencies
---------------
First install the following dependencies:
    [fftw3](http://fftw.org/)
    [tbb](https://github.com/intel/tbb)

Installation
--------------
First clone the repository with
```bash
git clone https://github.com/kthr/fluid-registration.git
```
Then create a directory "build" inside the cloned repository
```bash
cd fluid-registration
mkdir build
cd build
```

Finally build the binary
```bash
cmake ..
make
```
