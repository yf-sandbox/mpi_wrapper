MPI-Wrapper
====

MPI-Wrapper is a header only C++11 library for Message Passing Interface (MPI).

### Notes:
**<font color='Red'>This is under development!</font>**
Changes of specification are likely to happen.

## Features
- It is not necessary to include an additional library except for the MPI library. 
	(Google Test is optionally used for unit testing.)
- It is possible to code briefly.

### Why not use Boost::MPI?
Boost::MPI use other Boost libraries internally. 
Internal dependency is complicated in order to install selectively. 
On the other hand, Boost is too large to install all header files.

## Progress Table
|Function name|Status|
|--|--|
|Testing: | - |
|Point-to-Point Communication: | :heavy_check_mark: |
|Derived Datatypes: | :heavy_check_mark: |
|Communicators: | Partial |
|Groups: | - |
|FILE I/O: | - |
|Synchronization: | Partial |
|Collective Communication: | Partial |

<!-- https://gist.github.com/citrusui/07978f14b11adada364ff901e27c7f61 -->
#### Detailed progress
<details>
<summary> Collective Communication: </summary>
<p>

|Function name|Status|
|--|--|
| broadcast: |:heavy_check_mark:|
| reduction: |:heavy_check_mark:|
| allreduction: |:heavy_check_mark:|
| scatter: |-|
| gather: |-|
| alltoall: |-|
| allgather: |:heavy_check_mark:|

</p>
</details>

## Installation
#### Using CMake
```
mkdir build	# For building out of source
cd build
cmake .. 	# Generate native build script
cmake --build .	# Build test
ctest -V	# Execute testing
```
If you skip tests, you should replace the generation command of build script with:
```
cmake -DENABLE_TESTING=OFF ..
```

#### Without CMake
Just copy `mpi_wrapper.hpp` in your project and then build your project such as following:
```
mpic++ -std=c++11 your_source.cpp
```

## Usage
Data members of the derived type are indicated by serializer.
Therefore, if you use the derived type, it is necessary to define the serialize function.
(I hope that C/C++ support reflection.)
```
#include <mpi_wrapper.hpp>

/* The serialize function is defined inside a class. */
struct Test_Inside {
   int a;
   double b;
   template <class Archive> void serialize(Archive &archive) const {archive(a, b);}
   template <class Archive> void serialize(Archive &archive) {archive(a, b);}
};

/* The serialize function is defined outside a class. */
struct Test_Outside {
   int a;
   double b;
};
template <class Archive> void serialize(Archive &archive, const Test_Outside &t) {archive(t.a, t.b);}
template <class Archive> void serialize(Archive &archive, Test_Outside &t) {archive(t.a, t.b);}

int main(int argc, char *argv[]) {
   mpiw::initialize(&argc, &argv);
   mpiw::communicator world;
   
   const int size = world.size(); // Assume size=4
   const int rank = world.rank();
   
   /* Test_Inside and Test_Outside are compatible types. */
   struct Test_Inside ti = {1, 2.3};
   struct Test_Outside to;
   
   mpiw:allreduce(ti, to, mpiw::op::sum);
                  
   std::cout << "to.a: " << to.a << ", to.b: " << to.b << std::endl; // to.a: 4, to.b: 9.2

   return 0;
}
```

## Requirements
- OpenMPI (or MPICH or MVAPICH)
- C++11

## License
[MIT License](https://github.com/tcnksm/tool/blob/master/LICENCE)

