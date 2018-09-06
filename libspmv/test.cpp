#include <dlfcn.h>

#include <iostream>
#include <tuple>
#include <vector>

using harness_t = void(double *ov, double *a, double *iv, int *rowstr, int *colidx, int *rows);
using f_harness_t = void(float *ov, float *a, float *iv, int *rowstr, int *colidx, int *rows);

void usage(int argc, char **argv)
{
  if(argc != 2) {
    std::cerr << "Usage: " << argv[0] << " [version]\n";
    std::exit(1);
  }
}

std::string lib_path(const std::string& name)
{
  return "./" + name + ".so";
}

template <typename Floating, typename Harness>
bool test(Harness impl)
{
  auto A = std::vector<Floating>{5.0, 8.0, 3.0, 6.0};
  auto rowstr = std::vector<int>{1, 1, 3, 4, 5};
  auto colidx = std::vector<int>{1, 2, 3, 2};
  auto x = std::vector<Floating>{1.0, 2.0, 3.0, 4.0};
  auto y = std::vector<Floating>{0.0, 21.0, 9.0, 12.0};
  auto output = std::vector<Floating>(4);
  int rows = 4;
  
  impl(output.data(), A.data(), x.data(), rowstr.data(), colidx.data(), &rows);

  return std::equal(output.begin(), output.end(), y.begin());
}

int main(int argc, char **argv)
{
  usage(argc, argv);
  auto path = lib_path(argv[1]);

  void *lib = dlopen(path.c_str(), RTLD_NOW);
  if(!lib) {
    std::cerr << "Couldn't load " << path << '\n';
    std::cerr << dlerror() << '\n';
    std::exit(2);
  }

  dlerror();

  void *sym = dlsym(lib, "spmv_harness_");
  void *f_sym = dlsym(lib, "f_spmv_harness_");

  char *err = dlerror();
  if(err) {
    std::cerr << "Error loading SPMV implementation:\n" << err << '\n';
    std::exit(3);
  }

  auto harness = reinterpret_cast<harness_t *>(sym);
  auto f_harness = reinterpret_cast<f_harness_t *>(f_sym);

  //if(!test<double, harness_t>(harness) || !test<float, f_harness_t>(f_harness)) {
  if(!test<float, f_harness_t>(f_harness)) {
    std::cerr << "failure\n";
  } else {
    std::cerr << "success!\n";
  }
}
