#include <gtest/gtest.h>
#include <mpi_wrapper.hpp>
#include <iostream>
#include <csignal>
#include <cmath>

#include <random>
#include <algorithm>

// We have to write main function by myself for recieve the arguments in MPI application.

TEST(Test, Initialization) {
  EXPECT_EQ(1, 1);
  EXPECT_EQ(2, 2);
}

TEST(DISABLED_Test, Failure) {
  EXPECT_EQ(1, 1);
  EXPECT_EQ(2, 30000);
}

TEST(Test, Communicator) {
  mpiw::Communicator world;
  ASSERT_TRUE(0 <= world.rank());
  ASSERT_TRUE(4 == world.size()) << "size: " << world.size();
}

struct Inner {
  double v1;
  int v2;
  template <class Archive> void serialize(Archive &a) const {a(v1);a(v2);}
  template <class Archive> void serialize(Archive &a) {a(v1);a(v2);}
};
Inner operator+(const Inner &lhs, const Inner &rhs) {return Inner{lhs.v1+rhs.v1, lhs.v2+rhs.v2};}
struct Outer {
  char v1;
  Inner v2;
};
template <class Archive> void serialize(Archive &a, const Outer &o) {a(o.v1);a(o.v2);}
template <class Archive> void serialize(Archive &a, Outer &o) {a(o.v1);a(o.v2);}
Outer operator+(const Outer &lhs, const Outer &rhs) {
  return Outer{static_cast<char>(lhs.v1+rhs.v1), lhs.v2+rhs.v2};
}

TEST(Test, Derived_Datatype) {
  mpiw::Communicator world;
  const int rank = world.rank();
  const int size = world.size();
  Outer value = {static_cast<char>(rank+1), {0.2*(rank+1), 3*(rank+1)}};
  // Allreduce
  Outer sum = value;
  world.allreduce(value, sum, mpiw::op::plus);
  const int tri_num = size*(size+1)/2; // n-th triangular number (n=size)
  ASSERT_TRUE(tri_num == static_cast<int>(sum.v1)) << "sum.v1(rank"<<rank<<"): " << static_cast<int>(sum.v1);
  ASSERT_TRUE(0.2*tri_num == sum.v2.v1) << "sum.v2.v1(rank"<<rank<<"): " << sum.v2.v1;
  ASSERT_TRUE(3*tri_num == sum.v2.v2) << "sum.v2.v2(rank"<<rank<<"): " << sum.v2.v2;
}

struct Test1 {
  int a;
  double b;
  template <class Archive> void serialize(Archive &archive) {archive(a, b);}
  template <class Archive> void serialize(Archive &archive) const {archive(a, b);}
};
Test1 operator+(const Test1 &lhs, const Test1 &rhs) {return Test1{lhs.a+rhs.a, lhs.b+rhs.b};}
struct Test2 {
  int a;
  double b;
};
template <class Archive> void serialize(Archive &archive, Test2 &t) {archive(t.a, t.b);}
template <class Archive> void serialize(Archive &archive, const Test2 &t) {archive(t.a, t.b);}

TEST(Test, Compatible_Datatype) {
  mpiw::Communicator world;
  const int rank = world.rank();
  const int size = world.size();
  Test1 t1 = {1, 2.3};
  Test2 t2;
  world.allreduce(t1,t2,mpiw::op::plus);
  ASSERT_TRUE(t2.a == t1.a*size) << "t2.a(rank"<<rank<<"): " << t2.a;
  ASSERT_TRUE(t2.b == t1.b*size) << "t2.b(rank"<<rank<<"): " << t2.b;
}


TEST(Test, Send_Recv) {
  mpiw::Communicator world;
  const int rank = world.rank();
  const int size = world.size();
  ASSERT_TRUE(2 <= world.size()) << "size(rank"<<rank<<"): " << world.size();
  
  const int seed = 2019+rank; // Seed of rank0 is only use.
  int value = seed;

  // Broadcast
  const int depth = static_cast<int>(std::log(size)/std::log(2));
  for(int i=0; i<depth; ++i) {
    const int p2i = pow(2,i);
    const int target = (rank%(p2i*2)<p2i) ? rank+p2i : rank-p2i;
    if(target > rank) {
      world.send(value, target, 0);
    } else {
      world.recv(value, target, 0);
    }
  }
  // Varify
  if(rank == 0) {
    ASSERT_TRUE(value == seed) << "value(rank"<<rank<<"): " << value;
  } else {
    ASSERT_TRUE(value != seed) << "value(rank"<<rank<<"): " << value;
  }

  // Allreduce (using Butterfly reduction)
  int sum = value;
  for(int i=0; i<depth; ++i) {
    const int p2i = pow(2,i);
    const int target = (rank%(p2i*2)<p2i) ? rank+p2i : rank-p2i;
    int tmp = 0;
    if(target > rank) {
      world.recv(tmp, target, 0);
      world.send(sum, target, 0);
    } else {
      world.send(sum, target, 0);
      world.recv(tmp, target, 0);
    }
    sum += tmp;
  }
  // Varify
  if(rank == 0) {
    ASSERT_TRUE(sum == seed*size) << "sum(rank"<<rank<<"): " << sum;
  } else {
    ASSERT_TRUE(sum != seed*size) << "sum(rank"<<rank<<"): " << sum;
  }
}

TEST(Test, ISend_IRecv) {
  mpiw::Communicator world;
  const int rank = world.rank();
  const int size = world.size();
  ASSERT_TRUE(2 <= world.size()) << "size(rank"<<rank<<"): " << world.size();
  
  const int seed = 2019+rank; // Seed of rank0 is only use.
  int value = seed;

  // Broadcast
  const int depth = static_cast<int>(std::log(size)/std::log(2));
  for(int i=0; i<depth; ++i) {
    const int p2i = pow(2,i);
    const int target = (rank%(p2i*2)<p2i) ? rank+p2i : rank-p2i;
    if(target > rank) {
      world.isend(value, target, 0).wait();
    } else {
      world.irecv(value, target, 0).wait();
    }
  }
  // Varify
  if(rank == 0) {
    ASSERT_TRUE(value == seed) << "value(rank"<<rank<<"): " << value;
  } else {
    ASSERT_TRUE(value != seed) << "value(rank"<<rank<<"): " << value;
  }

  // Allreduce (using Butterfly reduction)
  int sum = value;
  for(int i=0; i<depth; ++i) {
    const int p2i = pow(2,i);
    const int target = (rank%(p2i*2)<p2i) ? rank+p2i : rank-p2i;
    int tmp = 0;
    mpiw::Request req_send = world.isend(sum, target, 0);
    mpiw::Request req_recv = world.irecv(tmp, target, 0);
    req_send.wait();
    req_recv.wait();
    sum += tmp;
  }
  // Varify
  if(rank == 0) {
    ASSERT_TRUE(sum == seed*size) << "sum(rank"<<rank<<"): " << sum;
  } else {
    ASSERT_TRUE(sum != seed*size) << "sum(rank"<<rank<<"): " << sum;
  }
}

TEST(Test, Broadcast_Allreduce) {
  mpiw::Communicator world;
  const int rank = world.rank();
  const int size = world.size();
  ASSERT_TRUE(2 <= world.size()) << "size(rank"<<rank<<"): " << world.size();
  
  const int seed = 2019+rank; // Seed of rank0 is only use.
  int value = seed;

  // Broadcast
  world.broadcast(value, 0);
  // Varify
  if(rank == 0) {
    ASSERT_TRUE(value == seed) << "value(rank"<<rank<<"): " << value;
  } else {
    ASSERT_TRUE(value != seed) << "value(rank"<<rank<<"): " << value;
  }

  // Allreduce
  int sum = value;
  world.allreduce(value, sum, mpiw::op::plus);
  // Varify
  if(rank == 0) {
    ASSERT_TRUE(sum == seed*size) << "sum(rank"<<rank<<"): " << sum;
  } else {
    ASSERT_TRUE(sum != seed*size) << "sum(rank"<<rank<<"): " << sum;
  }

  //  world.allreduce(seed, sum, mpiw::op::plus);
  //  std::cout<< "sum("<<rank<<"): " << sum << std::endl;
}

TEST(Test, IBroadcast_IAllreduce) {
  mpiw::Communicator world;
  const int rank = world.rank();
  const int size = world.size();
  ASSERT_TRUE(2 <= world.size()) << "size(rank"<<rank<<"): " << world.size();
  
  const int seed = 2019+rank; // Seed of rank0 is only use.
  int value = seed;

  // Broadcast
  world.ibroadcast(value, 0).wait();
  // Varify
  if(rank == 0) {
    ASSERT_TRUE(value == seed) << "value(rank"<<rank<<"): " << value;
  } else {
    ASSERT_TRUE(value != seed) << "value(rank"<<rank<<"): " << value;
  }

  // Allreduce
  int sum = value;
  world.iallreduce(value, sum, mpiw::op::plus).wait();
  // Varify
  if(rank == 0) {
    ASSERT_TRUE(sum == seed*size) << "sum(rank"<<rank<<"): " << sum;
  } else {
    ASSERT_TRUE(sum != seed*size) << "sum(rank"<<rank<<"): " << sum;
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  mpiw::initialize(&argc, &argv);
  return RUN_ALL_TESTS();
}
