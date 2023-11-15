#include "tp_image_utils_functions/GaussBlurEngine.h"
#include "tp_image_utils_functions/BoxBlur.h"


#include <vector>
#include <cmath>
#include <thread>
#include <atomic>
#include <sstream>

#ifdef USE_OPENCL
#include <CL/cl.hpp>

#define CaseReturnString(x) case x: return #x;

namespace
{

const char *opencl_errstr(cl_int err)
{
    switch (err)
    {
        CaseReturnString(CL_SUCCESS                        )
        CaseReturnString(CL_DEVICE_NOT_FOUND               )
        CaseReturnString(CL_DEVICE_NOT_AVAILABLE           )
        CaseReturnString(CL_COMPILER_NOT_AVAILABLE         )
        CaseReturnString(CL_MEM_OBJECT_ALLOCATION_FAILURE  )
        CaseReturnString(CL_OUT_OF_RESOURCES               )
        CaseReturnString(CL_OUT_OF_HOST_MEMORY             )
        CaseReturnString(CL_PROFILING_INFO_NOT_AVAILABLE   )
        CaseReturnString(CL_MEM_COPY_OVERLAP               )
        CaseReturnString(CL_IMAGE_FORMAT_MISMATCH          )
        CaseReturnString(CL_IMAGE_FORMAT_NOT_SUPPORTED     )
        CaseReturnString(CL_BUILD_PROGRAM_FAILURE          )
        CaseReturnString(CL_MAP_FAILURE                    )
        CaseReturnString(CL_MISALIGNED_SUB_BUFFER_OFFSET   )
        CaseReturnString(CL_COMPILE_PROGRAM_FAILURE        )
        CaseReturnString(CL_LINKER_NOT_AVAILABLE           )
        CaseReturnString(CL_LINK_PROGRAM_FAILURE           )
        CaseReturnString(CL_DEVICE_PARTITION_FAILED        )
        CaseReturnString(CL_KERNEL_ARG_INFO_NOT_AVAILABLE  )
        CaseReturnString(CL_INVALID_VALUE                  )
        CaseReturnString(CL_INVALID_DEVICE_TYPE            )
        CaseReturnString(CL_INVALID_PLATFORM               )
        CaseReturnString(CL_INVALID_DEVICE                 )
        CaseReturnString(CL_INVALID_CONTEXT                )
        CaseReturnString(CL_INVALID_QUEUE_PROPERTIES       )
        CaseReturnString(CL_INVALID_COMMAND_QUEUE          )
        CaseReturnString(CL_INVALID_HOST_PTR               )
        CaseReturnString(CL_INVALID_MEM_OBJECT             )
        CaseReturnString(CL_INVALID_IMAGE_FORMAT_DESCRIPTOR)
        CaseReturnString(CL_INVALID_IMAGE_SIZE             )
        CaseReturnString(CL_INVALID_SAMPLER                )
        CaseReturnString(CL_INVALID_BINARY                 )
        CaseReturnString(CL_INVALID_BUILD_OPTIONS          )
        CaseReturnString(CL_INVALID_PROGRAM                )
        CaseReturnString(CL_INVALID_PROGRAM_EXECUTABLE     )
        CaseReturnString(CL_INVALID_KERNEL_NAME            )
        CaseReturnString(CL_INVALID_KERNEL_DEFINITION      )
        CaseReturnString(CL_INVALID_KERNEL                 )
        CaseReturnString(CL_INVALID_ARG_INDEX              )
        CaseReturnString(CL_INVALID_ARG_VALUE              )
        CaseReturnString(CL_INVALID_ARG_SIZE               )
        CaseReturnString(CL_INVALID_KERNEL_ARGS            )
        CaseReturnString(CL_INVALID_WORK_DIMENSION         )
        CaseReturnString(CL_INVALID_WORK_GROUP_SIZE        )
        CaseReturnString(CL_INVALID_WORK_ITEM_SIZE         )
        CaseReturnString(CL_INVALID_GLOBAL_OFFSET          )
        CaseReturnString(CL_INVALID_EVENT_WAIT_LIST        )
        CaseReturnString(CL_INVALID_EVENT                  )
        CaseReturnString(CL_INVALID_OPERATION              )
        CaseReturnString(CL_INVALID_GL_OBJECT              )
        CaseReturnString(CL_INVALID_BUFFER_SIZE            )
        CaseReturnString(CL_INVALID_MIP_LEVEL              )
        CaseReturnString(CL_INVALID_GLOBAL_WORK_SIZE       )
        CaseReturnString(CL_INVALID_PROPERTY               )
        CaseReturnString(CL_INVALID_IMAGE_DESCRIPTOR       )
        CaseReturnString(CL_INVALID_COMPILER_OPTIONS       )
        CaseReturnString(CL_INVALID_LINKER_OPTIONS         )
        CaseReturnString(CL_INVALID_DEVICE_PARTITION_COUNT )
        default: return "Unknown OpenCL error code";
    }
}

const char kernel_boxBlurH_4[] =
R"V0G0N(

float3 getf3(float* v, int index){
  index *= 3;
  return (float3)(v[index], v[index+1], v[index+2]);
}

void setf3(float* dst, int index, float3 v){
  index *= 3;
  dst[index+0] = v.x;
  dst[index+1] = v.y;
  dst[index+2] = v.z;
}

void kernel boxBlurH_4(global float* scl_, global float* tcl_, int w, int h, int r)
{
  const int nCnl = 3;
  const int nHead = r+1;
  const int nTail = r;
  const int nBody = w-(nHead+nTail);

  const float iarr = 1.0f / (float)(r + r + 1);

  int const i = get_global_id(0);

  float* scl = scl_ + w * i* nCnl;
  float* tcl = tcl_ + w * i* nCnl;

  float* liSCL = scl;
  float* riSCL = scl+r*nCnl;
  float* tiTCL = tcl;

  float3 fv  = getf3(scl,0);
  float3 lv  = getf3(scl, w - 1);
  float3 val = (r + 1) * fv;

  for(size_t j=0; j<r; j++)
    val += getf3(scl, j);

  for(float* tiTCLMax=tiTCL+nHead*nCnl; tiTCL<tiTCLMax; riSCL+=nCnl, tiTCL+=nCnl)
  {
    val += getf3(riSCL,0) - fv;
    setf3(tiTCL, 0 , val*iarr);
  }


  for(float* tiTCLMax=tiTCL+nBody*nCnl; tiTCL<tiTCLMax; liSCL+=nCnl, riSCL+=nCnl, tiTCL+=nCnl)
  {
    val += getf3(riSCL,0) - getf3(liSCL,0);
    setf3(tiTCL, 0 , val*iarr);
  }


  for(float* tiTCLMax=tiTCL+nTail*nCnl; tiTCL<tiTCLMax; liSCL+=nCnl, tiTCL+=nCnl)
  {
    val += lv - getf3(liSCL,0);
    setf3(tiTCL, 0 , val*iarr);
  }
}
)V0G0N";



const char  kernel_boxBlurT_4[] =
R"V0G0N(
void kernel boxBlurT_4(global float* scl_, global float* tcl_, int w, int h, int r)
{
  const int nCnl = 3;
  const float iarr = 1.0f / (r + r + 1);

  size_t const i = get_global_id(0);

  float* scl = scl_ + i * nCnl;
  float* tcl = tcl_ + i * nCnl;

  size_t ti = 0;
  size_t li = 0;
  size_t ri = r * w;

  float3 fv = getf3(scl, 0);
  float3 lv = getf3(scl, w * (h - 1));
  float3 val = (r + 1)*fv;

  for(size_t j=0; j<r; j++)
    val += getf3(scl, j * w);

  for(size_t j=0; j<=r; j++)
  {
    val += getf3(scl, ri) - fv;
    setf3(tcl, ti,val*iarr);
    ri += w; ti += w;
  }

  for(size_t j=r+1; j<h-r; j++)
  {
    val += getf3(scl, ri) - getf3(scl, li);
    setf3(tcl,ti, val*iarr);
    li += w;
    ri += w;
    ti += w;
  }

  for(size_t j=h-r; j<h; j++)
  {
    val += lv-getf3(scl,li);
    setf3(tcl,ti, val*iarr);
    li += w;
    ti += w;
  }
}
)V0G0N";


} //namespace
#endif


namespace tp_image_utils_functions
{

#ifdef USE_OPENCL

class GaussBlurEngine::GaussBlurAccelerator
{

private:
  std::vector<cl::Platform> all_platforms;
  std::vector<cl::Device> all_devices;
  cl::Device default_device;
  cl::Context context;
  cl::Program program;
  cl::Kernel boxBlurH_4;
  cl::Kernel boxBlurT_4;
  cl::CommandQueue queue;

  std::stringstream errors;
  std::stringstream info;

public:

  GaussBlurAccelerator()
  {
    cl::Platform::get(&all_platforms);
    if(all_platforms.size()==0){
      errors << "No platforms found. Check OpenCL installation!" << std::endl;
      return;
    }
    cl::Platform default_platform=all_platforms[0];
    info << "Using platform: "<< default_platform.getInfo<CL_PLATFORM_NAME>()<< std::endl;

    //get default device of the default platform
    default_platform.getDevices(CL_DEVICE_TYPE_ALL, &all_devices);
    if(all_devices.size()==0){
      errors << " No devices found. Check OpenCL installation!" << std::endl;
      return;
    }

    default_device=all_devices[0];
    info << "Using device: "<< default_device.getInfo<CL_DEVICE_NAME>()<< std::endl;

    context = cl::Context({default_device});

    //Next we need to create the program which we want to execute on our device:
    cl::Program::Sources sources;

    //Next we need our kernel sources to build. We also check for the errors at building:
    sources.push_back({kernel_boxBlurH_4, sizeof(kernel_boxBlurH_4)-1});
    sources.push_back({kernel_boxBlurT_4, sizeof(kernel_boxBlurT_4)-1});

    program = cl::Program(context, sources);
    if(program.build({default_device})!=CL_SUCCESS){
      errors << " Error building: "<<program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(default_device)<< std::endl;
      return;
    }

    boxBlurH_4 = cl::Kernel(program, "boxBlurH_4");
    boxBlurT_4 = cl::Kernel(program, "boxBlurT_4");

    //create queue to which we will push commands for the device.
    queue = cl::CommandQueue(context,default_device);
  }

  void doBlur(float* scl, size_t w, size_t h, size_t r)
  {
    if(!isFailed())
    {
      //For arrays source, target, we need to allocate the space on the device:
      size_t bufSize = sizeof(float)*h*w*3;

      // create buffers on the device
      cl::Buffer buffer_A(context,CL_MEM_READ_WRITE, bufSize);
      cl::Buffer buffer_B(context,CL_MEM_READ_WRITE, bufSize);

      //write arrays A and B to the device
      queue.enqueueWriteBuffer(buffer_A,CL_TRUE,0, bufSize, scl);
      //queue.enqueueWriteBuffer(buffer_B,CL_TRUE,0, bufSize, tcl);
      queue.finish();

      std::vector<size_t> bxs = tp_image_utils_functions::boxesForGauss(float(r), 3);

      for(int i:{0,1,2})
      {
        if(1)
        {
          boxBlurH_4.setArg(0, buffer_A);
          boxBlurH_4.setArg(1, buffer_B);
          boxBlurH_4.setArg(2, int(w));
          boxBlurH_4.setArg(3, int(h));
          boxBlurH_4.setArg(4, int(((bxs[i] - 1) / 2)));
          queue.enqueueNDRangeKernel(boxBlurH_4, cl::NullRange, cl::NDRange(h), cl::NullRange);
          auto err_code = queue.finish();
          if(err_code !=CL_SUCCESS)
          {
            errors <<" Error executing: " << opencl_errstr(err_code) << std::endl;
            return;
          }
        }

        if(1)
        {
          boxBlurT_4.setArg(0, buffer_B);
          boxBlurT_4.setArg(1, buffer_A);
          boxBlurT_4.setArg(2, int(w));
          boxBlurT_4.setArg(3, int(h));
          boxBlurT_4.setArg(4, int(((bxs[i] - 1) / 2)));
          queue.enqueueNDRangeKernel(boxBlurT_4, cl::NullRange, cl::NDRange(w), cl::NullRange);
          auto err_code = queue.finish();
          if(err_code !=CL_SUCCESS){
            errors <<" Error executing: " << opencl_errstr(err_code) << std::endl;
            return;
          }
        }
        //std::swap(buffer_B,buffer_A);
      }

      // please note that result is in array_A as we alternalte input and output arrays
      queue.enqueueReadBuffer(buffer_A, CL_TRUE, 0, bufSize, scl);
      queue.finish();
    }
  }

  bool isFailed()
  {
    return errors.tellp() != 0;
  }

  std::string getErrorString()
  {
    return errors.str();
  }

  std::string getInfoString()
  {
    return info.str();
  }
};
#else
class GaussBlurEngine::GaussBlurAccelerator{};
#endif

GaussBlurEngine::GaussBlurEngine(){}
GaussBlurEngine::~GaussBlurEngine(){}

void GaussBlurEngine::doBlur(float* scl, size_t w, size_t h, size_t r)
{
#ifdef USE_OPENCL
  if(!oclGaussBlur)
    oclGaussBlur.reset(new GaussBlurAccelerator());

  if(!oclGaussBlur->isFailed())
    oclGaussBlur->doBlur(scl, w, h, r);
  else
#endif
  {
    std::unique_ptr<float[]> aux(new float[w*h*3]);
    gaussBlur_4_cpu(scl, aux.get(), w, h, r);
  }
}

std::string GaussBlurEngine::getErrorString()
{
#ifdef USE_OPENCL
  if(oclGaussBlur && oclGaussBlur->isFailed())
    return oclGaussBlur->getErrorString();
  return "OpenCL capability not found using CPU for blur";
#else
  return "BUILT WITHOUT USE_OPENCL FLAG";
#endif
}

std::string GaussBlurEngine::getInfoString()
{
#ifdef USE_OPENCL
  if(oclGaussBlur)
    return oclGaussBlur->getInfoString();
  return "OpenCL capability not found using CPU for blur";
#else
  return "BUILT WITHOUT USE_OPENCL FLAG";
#endif
}

}
