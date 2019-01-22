#define CL_HPP_ENABLE_EXCEPTION
#define CL_HPP_MINIMUM_OPENCL_VERSION BUILD_CLVERSION
#define CL_HPP_TARGET_OPENCL_VERSION BUILD_CLVERSION

#include <CL/cl.hpp>

#include <iostream>
#include <vector>

int main()
{
  auto platforms = std::vector<cl::Platform>{};
  auto platforms_stat = cl::Platform::get(&platforms);
  if(platforms_stat != CL_SUCCESS) {
    std::cerr << "Couldn't get OpenCL platforms"
              << " [" << platforms_stat << "]\n";
    return 1;
  }

  for(auto i = 0; i < platforms.size(); ++i) {
    auto const& p = platforms.at(i);

    auto name = p.getInfo<CL_PLATFORM_NAME>();
    auto vendor = p.getInfo<CL_PLATFORM_VENDOR>();

    std::cout << "Platform " << i << ": " << name << "\n";

    auto devices = std::vector<cl::Device>{};
    auto devices_stat = p.getDevices(CL_DEVICE_TYPE_ALL, &devices);
    if(devices_stat != CL_SUCCESS) {
      std::cerr << "Couldn't get OpenCL devices"
                << " [" << devices_stat << "]\n";
    }

    for(auto i = 0; i < devices.size(); ++i) {
      auto const& dev = devices.at(i);

      auto name = dev.getInfo<CL_DEVICE_NAME>();
      auto avail = dev.getInfo<CL_DEVICE_AVAILABLE>() ? "available" : "unavailable";
      auto type = [&] {
        switch(dev.getInfo<CL_DEVICE_TYPE>()) {
          case CL_DEVICE_TYPE_CPU:
            return "CPU";
          case CL_DEVICE_TYPE_GPU:
            return "GPU";
          case CL_DEVICE_TYPE_ACCELERATOR:
            return "accelerator";
          default:
            return "?";
        }
      }();

      std::cout << "  Device " << i << ": " 
                << name << " (" << type << ", "
                << avail << ")\n";
    }
  }
}
