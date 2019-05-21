#pragma once

#include <mpi.h>
#include <cstdlib>

#include <iostream>

namespace mpiw {
  void initialize(int *ptr_argc = nullptr, char **ptr_argv[] = nullptr);

  class Archive {
  public:
    Archive(MPI_Datatype *type);
    ~Archive();
    void parse() {/* Do nothing */}
    template <class Type, class... Types>
    void parse(Type && type, Types && ...args);
    template <class... Types>
    Archive &operator()(Types && ...args);
    template <class Type> void head(const Type &v);
    template <class T> void eval(const T &v);
    template <class T, int N> void eval(const T (&v)[N]);
  private:
    MPI_Datatype *new_type;
    const char *head_ptr;
    std::vector<int> lengths;
    std::vector<MPI_Aint> offsets;
    std::vector<MPI_Datatype> types;
    Archive() = delete;
  };

  class CustomOperator {
  public:
    CustomOperator(MPI_User_function *func_ptr):op() {
      MPI_Op_create(func_ptr, 1, &op);
    }
    ~CustomOperator() {MPI_Op_free(&op);}
    operator MPI_Op() const {return op;}
  private:
    MPI_Op op;
    CustomOperator() = delete;
  };

  /* swich for sfinae */
  template <class T> using Enabler_t = typename std::enable_if<T::value,std::nullptr_t>::type;
  template <class T> using Disabler_t = typename std::enable_if<!T::value,std::nullptr_t>::type;

  namespace op {
    namespace factory {
      struct maximum {
	template <class T, Enabler_t<std::is_arithmetic<T>> =nullptr> MPI_Op generate() const {
	  return MPI_MAX;
	}
	template <class T, Disabler_t<std::is_arithmetic<T>> =nullptr> MPI_Op generate() const {
	  static const CustomOperator op(reinterpret_cast<MPI_User_function *>(maximum::impl<T>));
	  return op;
	}
	template <class T> static void impl(T *in, T *inout, int *len, MPI_Datatype *) {
	  for(int i=0; i<*len; ++i) { 
	    if(inout[i] < in[i]) { inout[i] = in[i]; }
	  }
	  return;
	}
      };
      struct minimum {
	template <class T, Enabler_t<std::is_arithmetic<T>> = nullptr> MPI_Op generate() const {
	  return MPI_MIN;
	}
	template <class T, Disabler_t<std::is_arithmetic<T>> = nullptr> MPI_Op generate() const {
	  static const CustomOperator op(reinterpret_cast<MPI_User_function *>(minimum::impl<T>));
	  return op;
	}
	template <class T> static void impl(T *in, T *inout, int *len, MPI_Datatype *) {
	  for(int i=0; i<*len; ++i) { 
	    if(inout[i] > in[i]) { inout[i] = in[i]; }
	  }
	  return;
	}
      };
      struct plus {
	template <class T, Enabler_t<std::is_arithmetic<T>> = nullptr> MPI_Op generate() const {
	  return MPI_SUM;
	}
	template <class T, Disabler_t<std::is_arithmetic<T>> = nullptr> MPI_Op generate() const {
	  static const CustomOperator op(reinterpret_cast<MPI_User_function *>(plus::impl<T>));
	  return op;
	}
	template <class T> static void impl(T *in, T *inout, int *len, MPI_Datatype *) {
	  for(int i=0; i<*len; ++i) {
	    //	  inout[i] += in[i];
	    inout[i] = inout[i] + in[i]; // fix using SFINAE
	  }
	  return;
	}
      };
      struct multiplies {
	template <class T, Enabler_t<std::is_arithmetic<T>> = nullptr> MPI_Op generate() const {
	  return MPI_PROD;
	}
	template <class T, Disabler_t<std::is_arithmetic<T>> = nullptr> MPI_Op generate() const {
	  static const CustomOperator op(reinterpret_cast<MPI_User_function *>(multiplies::impl<T>));
	  return op;
	}
	template <class T> static void impl(T *in, T *inout, int *len, MPI_Datatype *) {
	  for(int i=0; i<*len; ++i) {
	    //	  inout[i] *= in[i];
	    inout[i] = inout[i] * in[i]; // fix using SFINAE
	  }
	  return;
	}
      };
      struct logical_and {
	template <class T, Enabler_t<std::is_arithmetic<T>> = nullptr> MPI_Op generate() const {
	  return MPI_LAND;
	}
	template <class T, Disabler_t<std::is_arithmetic<T>> = nullptr> MPI_Op generate() const {
	  static const CustomOperator op(reinterpret_cast<MPI_User_function *>(logical_and::impl<T>));
	  return op;
	}
	template <class T> static void impl(T *in, T *inout, int *len, MPI_Datatype *) {
	  for(int i=0; i<*len; ++i) { inout[i] = inout[i] && in[i]; }
	  return;
	}
      };
      struct logical_or {
	template <class T, Enabler_t<std::is_arithmetic<T>> = nullptr> MPI_Op generate() const {
	  return MPI_LOR;
	}
	template <class T, Disabler_t<std::is_arithmetic<T>> = nullptr> MPI_Op generate() const {
	  static const CustomOperator op(reinterpret_cast<MPI_User_function *>(logical_or::impl<T>));
	  return op;
	}
	template <class T> static void impl(T *in, T *inout, int *len, MPI_Datatype *) {
	  for(int i=0; i<*len; ++i) { inout[i] = inout[i] || in[i]; }
	  return;
	}
      };
      struct logical_xor {
	template <class T, Enabler_t<std::is_arithmetic<T>> = nullptr> MPI_Op generate() const {
	  return MPI_LXOR;
	}
	template <class T, Disabler_t<std::is_arithmetic<T>> = nullptr> MPI_Op generate() const {
	  static const CustomOperator op(reinterpret_cast<MPI_User_function *>(logical_xor::impl<T>));
	  return op;
	}
	template <class T> static void impl(T *in, T *inout, int *len, MPI_Datatype *) {
	  for(int i=0; i<*len; ++i) { inout[i] = !(inout[i]) != !(in[i]); }
	  return;
	}
      };
      struct bitwise_and {
	template <class T, Enabler_t<std::is_arithmetic<T>> = nullptr> MPI_Op generate() const {
	  return MPI_BAND;
	}
	template <class T, Disabler_t<std::is_arithmetic<T>> = nullptr> MPI_Op generate() const {
	  static const CustomOperator op(reinterpret_cast<MPI_User_function *>(bitwise_and::impl<T>));
	  return op;
	}
	template <class T> static void impl(T *in, T *inout, int *len, MPI_Datatype *) {
	  for(int i=0; i<*len; ++i) {
	    //	    inout[i] &= in[i];
	    inout[i] = inout[i] & in[i]; // fix using SFINAE
	  }
	  return;
	}
      };
      struct bitwise_or {
	template <class T, Enabler_t<std::is_arithmetic<T>> = nullptr> MPI_Op generate() const {
	  return MPI_BOR;
	}
	template <class T, Disabler_t<std::is_arithmetic<T>> = nullptr> MPI_Op generate() const {
	  static const CustomOperator op(reinterpret_cast<MPI_User_function *>(bitwise_or::impl<T>));
	  return op;
	}
	template <class T> static void impl(T *in, T *inout, int *len, MPI_Datatype *) {
	  for(int i=0; i<*len; ++i) {
	    //	    inout[i] |= in[i];
	    inout[i] = inout[i] | in[i]; // fix using SFINAE
	  }
	  return;
	}
      };
      struct bitwise_xor {
	template <class T, Enabler_t<std::is_arithmetic<T>> = nullptr> MPI_Op generate() const {
	  return MPI_BXOR;
	}
	template <class T, Disabler_t<std::is_arithmetic<T>> = nullptr> MPI_Op generate() const {
	  static const CustomOperator op(reinterpret_cast<MPI_User_function *>(bitwise_xor::impl<T>));
	  return op;
	}
	template <class T> static void impl(T *in, T *inout, int *len, MPI_Datatype *) {
	  for(int i=0; i<*len; ++i) {
	    //	    inout[i] ^= in[i];
	    inout[i] = inout[i] ^ in[i]; // fix using SFINAE
	  }
	  return;
	}
      };
    }
    const factory::maximum maximum;
    const factory::minimum minimum;
    const factory::plus plus;
    const factory::multiplies multiplies;
    const factory::logical_and logical_and;
    const factory::logical_or logical_or;
    const factory::logical_xor logical_xor;
    const factory::bitwise_and bitwise_and;
    const factory::bitwise_or bitwise_or;
    const factory::bitwise_xor bitwise_xor;
  }


  
  template <typename T> MPI_Datatype mpi_type(const T &v);

  class CustomDatatype {
  public:
    template <typename T> CustomDatatype(const T &v):type() {
      { // local scope to call destractor of Archive class
	Archive arch(&type);
	arch.head(v);
	serialize_handler(arch,v);
      }
      MPI_Type_commit(&type);
    }
    ~CustomDatatype() {MPI_Type_free(&type);}
    operator MPI_Datatype() const {return type;}
  private:
    MPI_Datatype type;
    CustomDatatype() = delete;
  };

  template <class T>
  struct DataInfo {
    DataInfo(T &v):ptr(&v),count(1),type(mpi_type(v)) {}
    template <int N> DataInfo(T (&v)[N]):ptr(&v[0]),count(N),type(mpi_type(v[0])) {}
    template <template <class,class=std::allocator<T>> class C> DataInfo(C<T> &v):ptr(v.data()),count(v.size()),type(mpi_type(v[0])) {}
    ~DataInfo() {}
    DataInfo() = delete;
    void *ptr;
    int count;
    MPI_Datatype type;
  };

  template <class T>
  struct ConstDataInfo {
    ConstDataInfo(const T &v):ptr(&v),count(1),type(mpi_type(v)) {}
    ConstDataInfo(T &&v):ptr(&v),count(1),type(mpi_type(v)) {}
    template <int N> ConstDataInfo(const T (&v)[N]):ptr(&v[0]),count(N),type(mpi_type(v[0])) {}
    template <template <class, class=std::allocator<T>> class C> ConstDataInfo(const C<T> &v):ptr(v.data()),count(v.size()),type(mpi_type(v[0])) {}
    ~ConstDataInfo() {}
    ConstDataInfo() = delete;
    const void *ptr;
    int count;
    MPI_Datatype type;
  };

  template <typename T> MPI_Datatype mpi_type(const T &v) {
    /* generate new datatype because of undefined datatype */
    static CustomDatatype t(v);
    return t;//.value();
  }
  /* explicitly specialize for primitive datatype */
  template <> MPI_Datatype mpi_type(const char &) { return MPI_CHAR; }
  template <> MPI_Datatype mpi_type(const signed char &) { return MPI_SIGNED_CHAR; }
  template <> MPI_Datatype mpi_type(const signed short int &) { return MPI_SHORT; }
  template <> MPI_Datatype mpi_type(const signed int &) { return MPI_INT; }
  template <> MPI_Datatype mpi_type(const signed long int &) { return MPI_LONG; }
  template <> MPI_Datatype mpi_type(const unsigned char &) { return MPI_UNSIGNED_CHAR; }
  template <> MPI_Datatype mpi_type(const unsigned short int &) { return MPI_UNSIGNED_SHORT; }
  template <> MPI_Datatype mpi_type(const unsigned int &) { return MPI_UNSIGNED; }
  template <> MPI_Datatype mpi_type(const unsigned long int &) { return MPI_UNSIGNED_LONG; }
  template <> MPI_Datatype mpi_type(const float &) { return MPI_FLOAT; }
  template <> MPI_Datatype mpi_type(const double &) { return MPI_DOUBLE; }
  template <> MPI_Datatype mpi_type(const long double &) { return MPI_LONG_DOUBLE; }

  template <typename T, class A> class has_serialize {
  public:
    static const bool value;
  private:
    template <typename U, typename V> static auto check(U u, V v) -> decltype(u.serialize(v), std::true_type());
    static auto check(...) -> std::false_type;
  };
  template <typename T, typename A> const bool has_serialize<T,A>::value = decltype(has_serialize<T,A>::check(std::declval<T>(), std::declval<A>()))::value;

  template <class Archive, typename T, Enabler_t<has_serialize<T,Archive>> = nullptr>
  auto serialize_handler(Archive &archive, T &&arg)->void {
    std::cout << "Has serialize()." << std::endl;
    return arg.serialize(archive);
  }
  template <class Archive, typename T, Disabler_t<has_serialize<T,Archive>> = nullptr>
  auto serialize_handler(Archive &archive, T &&arg)-> void {
    std::cout << "Dosen't have serialize()." << std::endl;
    return serialize(archive, arg);
  }

  Archive::Archive(MPI_Datatype *type):new_type(type),head_ptr(nullptr) {}
  Archive::~Archive() {MPI_Type_create_struct(lengths.size(), lengths.data(), offsets.data(), types.data(), new_type);}
  template <class Type, class... Types> void Archive::parse(Type && type, Types && ...args) {
    eval(type);
    parse(std::forward<Types>(args)...);
  }
  template <class... Types> Archive &Archive::operator()(Types && ...args) {
    parse(std::forward<Types>(args)...);
    return *this;
  }
  template <class Type> void Archive::head(const Type &v) {
    head_ptr = reinterpret_cast<const char*>(&v);
  }
  template <class T> void Archive::eval(const T &v) {
    lengths.push_back(1);
    offsets.push_back(reinterpret_cast<const char *>(&v)-head_ptr);
    types.push_back(mpi_type(v));
  }
  template <class T, int N> void Archive::eval(const T (&v)[N]) {
    lengths.push_back(N);
    offsets.push_back(reinterpret_cast<const char *>(v)-head_ptr);
    types.push_back(mpi_type(v[0])); // ToDo fix
  }

  class Status {
  public:
    Status(const MPI_Status &sta):sta(sta) {std::cout<<"status constract"<<std::endl;}
    ~Status() {std::cout<<"status constract"<<std::endl;}
    MPI_Status sta; // ToDo : move to private
  private:
    Status() = delete;
  };

  class Request {
  public:
    Request(const MPI_Request &req):req(req),sta() {/*std::cout << "request constract" << std::endl;*/}
    ~Request() {/*std::cout << "request destract" << std::endl;*/}
    /*** Status::operator bool() or not ***/
    bool wait() {
      if(MPI_Wait(&req, &sta) != MPI_SUCCESS) {
	exit(EXIT_FAILURE);
      }
      return true;
    }
    bool test() {
      int flag = 0;
      if(MPI_Test(&req, &flag, &sta) != MPI_SUCCESS) {
	exit(EXIT_FAILURE);
      }
      return (flag) ? true : false;
    }
    Status status() {
      return Status(sta);
    }
  private:
    MPI_Request req;
    MPI_Status sta;
    Request() = delete;
  };

  template <bool flag, class IsTrue, class IsFalse>
  struct choose;
  template <class IsTrue, class IsFalse>
  struct choose<true, IsTrue, IsFalse> {
    typedef IsTrue type;
  };
  template <class IsTrue, class IsFalse>
  struct choose<false, IsTrue, IsFalse> {
    typedef IsFalse type;
  };

  class Communicator {
  public:
    Communicator():comm(MPI_COMM_WORLD) {initialize();}
    int size() const {
      int x;
      if(MPI_Comm_size(comm, &x) != MPI_SUCCESS) {
	exit(EXIT_FAILURE);
      }
      return x;
    }
    int rank() const {
      int x;
      if(MPI_Comm_rank(comm, &x) != MPI_SUCCESS) {
	exit(EXIT_FAILURE);
      }
      return x;
    }

    /*** https://stackoverflow.com/questions/51609496/implicit-type-conversion-with-templated-function-parameters ***/
    template <class T, class... Args> void send(T &&t, Args&&... args) {
      send<typename std::decay<T>::type>(std::forward<T>(t), std::forward<Args>(args)...);
    }
    template <class T> void send(ConstDataInfo<T> &&data, const int &target, const int &tag) {
      MPI_Send(data.ptr, data.count, data.type, target, tag, comm);
      return;
    }
    template <class T, class... Args> void recv(T &&t, Args&&... args) {
      recv<typename std::decay<T>::type>(std::forward<T>(t), std::forward<Args>(args)...);
    }
    template <class T> void recv(DataInfo<T> &&data, const int &target, const int &tag) {
      MPI_Recv(data.ptr, data.count, data.type, target, tag, comm, MPI_STATUS_IGNORE);
      return;
    }
  
    template <class T, class... Args> Request isend(T &&t, Args&&... args) {
      return isend<typename std::decay<T>::type>(std::forward<T>(t), std::forward<Args>(args)...);
    }
    template <class T> Request isend(const ConstDataInfo<T> data, const int &target, const int &tag) {
      MPI_Request req;
      if(MPI_Isend(data.ptr, data.count, data.type, target, tag, comm, &req) != MPI_SUCCESS) {
	exit(EXIT_FAILURE);
      }
      return Request(req); // RVO
    }
    template <class T, class... Args> Request irecv(T &&t, Args&&... args) {
      return irecv<typename std::decay<T>::type>(std::forward<T>(t), std::forward<Args>(args)...);
    }
    template <class T> Request irecv(DataInfo<T> &&data, const int &target, const int &tag) {
      MPI_Request req;
      if(MPI_Irecv(data.ptr, data.count, data.type, target, tag, comm, &req) != MPI_SUCCESS) {
	exit(EXIT_FAILURE);
      }
      return Request(req); // RVO
    }
    template <class T> void broadcast(T &&t, const int &target) {
      return broadcast_impl<typename std::decay<T>::type>(std::forward<T>(t), target);
    }
    template <class T> void broadcast_impl(DataInfo<T> data, const int &target) {
      MPI_Bcast(data.ptr, data.count, data.type, target, comm);
      return;
    }
    template <class T, class U, class OpFactory> void allreduce(T &&t1, U &&t2, const OpFactory &of=op::plus) {
      return allreduce_impl<typename std::decay<T>::type, typename std::decay<U>::type>(std::forward<T>(t1), std::forward<U>(t2), of);
    }
    template <class T, class U, class OpFactory> void allreduce_impl(const ConstDataInfo<T> sdata, DataInfo<U> rdata, const OpFactory &of) {
      if(sdata.count != rdata.count) {
	std::cerr << "Warning(Allreduce): Data size is different between send and receive data." << std::endl;
      }
      if(sdata.type != rdata.type) {
	std::cerr << "Warning(Allreduce): Data type is different between send and receive data." << std::endl;
      }
      MPI_Allreduce(sdata.ptr, rdata.ptr, rdata.count, rdata.type, of.template generate<T>(), comm);
      return;
    }

    template <class T> Request ibroadcast(T &&t, const int &target) {
      return ibroadcast_impl<typename std::decay<T>::type>(std::forward<T>(t), target);
    }
    template <class T> Request ibroadcast_impl(DataInfo<T> &&data, const int &target) {
      MPI_Request req;
      MPI_Ibcast(data.ptr, data.count, data.type, target, comm, &req);
      return Request(req); // RVO
    }
    template <class T, class U, class OpFactory> Request iallreduce(T &&t, U &&t2, const OpFactory &of=op::plus) {
      return iallreduce_impl<typename std::decay<T>::type, typename std::decay<U>::type>(std::forward<T>(t), std::forward<U>(t2), of);
    }
    template <class T, class U, class OpFactory> Request iallreduce_impl(const ConstDataInfo<T> &sdata, DataInfo<T> &&rdata, const OpFactory &of) {
      MPI_Request req;
      if(sdata.count != rdata.count) {
	std::cerr << "Warning(IAllreduce): Data size is different between send and receive data." << std::endl;
      }
      if(sdata.type != rdata.type) {
	std::cerr << "Warning(IAllreduce): Data type is different between send and receive data." << std::endl;
      }
      MPI_Iallreduce(sdata.ptr, rdata.ptr, rdata.count, rdata.type, of.template generate<T>(), comm, &req);
      return Request(req); // RVO
    }

    void barrier() {
      MPI_Barrier(comm);
      return;
    }
    double wtime() {
      return MPI_Wtime();
    }
  private:
    MPI_Comm comm;
  };

  namespace {
    // implementation
    class mpi_impl {
      mpi_impl() {}
      ~mpi_impl() {}
    };

    class Environment {
    public:
      Environment(int *argc, char **argv[]) {MPI_Init(argc, argv);}
      ~Environment() {MPI_Finalize();}
    private:
      Environment() = delete;
    };
  }

  void initialize(int *ptr_argc, char **ptr_argv[]) {
    static Environment env(ptr_argc, ptr_argv);
    //  volatile static Environment env(ptr_argc, ptr_argv);
    //  (void) env;
  }

} // end namespace
