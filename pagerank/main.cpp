#include <mm/mm.h>

#include <dlfcn.h>

#include <algorithm>
#include <chrono>
#include <iostream>
#include <random>
#include <string>

class dynamic_library {
public:
  dynamic_library(const std::string& path)
  {
    lib_ = dlopen(path.c_str(), RTLD_NOW);
    if(!lib_) {
      throw std::runtime_error(dlerror());
    }
  }

  ~dynamic_library()
  {
    if(lib_) {
      dlclose(lib_);
    }
  }

  template <typename Func>
  Func *symbol(const std::string& symbol)
  {
    dlerror();
    void *sym = dlsym(lib_, symbol.c_str());
    char *err = dlerror();
    if(err) {
      throw std::runtime_error(err);
    }

    return reinterpret_cast<Func *>(sym);
  }

private:
  void *lib_;
};

struct benchmark_stats {
  std::string path;
  double time;
  size_t rows;
  size_t cols;
  size_t nnz;
  size_t iters;
};

std::ostream& operator<<(std::ostream& os, benchmark_stats const& stats)
{
  os << stats.path << " "
     << stats.time << " "
     << stats.rows << " "
     << stats.cols << " "
     << stats.nnz << " "
     << stats.iters;

  return os;
}

std::vector<double> random_starting_vector(size_t rows)
{
  std::vector<double> x(rows, 0);
  auto rd = std::random_device{};
  auto dist = std::uniform_real_distribution<>(0, 1.0);

  std::generate(x.begin(), x.end(), [&rd,&dist] {
    return dist(rd);
  });

  auto sum = std::accumulate(x.begin(), x.end(), 0.0);
  std::for_each(x.begin(), x.end(), [sum](auto& val) {
    val /= sum;
  });

  return x;
}

template <typename Func>
benchmark_stats run_benchmark(Func&& harness, std::string const& path)
{
  auto mat = mm::coordinate_matrix::read_from_file(path);
  if(mat.rows() != mat.cols()) {
    throw std::runtime_error("Pagerank needs a square matrix!");
  }
  mat.normalise();

  auto csr = mm::csr_matrix(mm::one_based_index, mat);
  auto d = 0.85;
  csr.scale(d);

  auto x = random_starting_vector(csr.rows());
  auto y = std::vector<double>(csr.rows(), 0);

  auto raw = mm::csr(csr);

  auto start = std::chrono::system_clock::now();

  volatile double error;
  std::vector<double> last_vector;

  auto iters = 1024u;
  for(auto i = 0; i < iters; ++i) {
    error = 0.0;

    last_vector = x;

    auto mean = std::accumulate(x.begin(), x.end(), 0.0) / x.size();
    auto add_term = (1 - d) * mean;

    /* v = (d * M) * v + ((1 - d)) * v.average(); */
    /* error = (v - last_v).l2norm(); */
    harness(y.data(), raw.a, x.data(), raw.rowstr, raw.colidx, raw.rows);

    std::for_each(y.begin(), y.end(), [add_term] (auto& out) {
      out += add_term;
    });

    x = y;
    for(auto i = 0; i < csr.rows(); ++i) {
      error += std::pow(x.at(i) - last_vector.at(i), 2);
    }
    error = std::sqrt(error);
  }

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> time = end - start;

  return { path, time.count(), csr.rows(), csr.cols(), csr.nnz(), iters };
}

using harness_t = void (double *, const double *, const double *, const int *, const int *, const int *);

void usage(char **argv)
{
  std::cerr << "Usage: " << argv[0] << " library" << " label" << " [benchmarks...]\n";
  std::cerr << "Arguments:\n";
  std::cerr << "  library: Path to a shared library with a compatible spmv_harness_ implementation\n";
  std::cerr << "  label: Text label to be printed at the end of each row to identify the implementation used\n";
  std::cerr << "  benchmarks: One or more paths to matrix market files to benchmark\n";
}

int main(int argc, char **argv)
{
  if(argc < 3) {
    usage(argv);
    return 1;
  }

  auto dyn = dynamic_library(argv[1]);
  auto harness = dyn.symbol<harness_t>("spmv_harness_");

  for(int i = 3; i < argc; ++i) {
    auto path = std::string(argv[i]);
    auto stats = run_benchmark(harness, path);
    std::cout << stats << " " << argv[2] << '\n';
  }
}
