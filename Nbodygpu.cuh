#include <cuda_runtime_api.h>
#include <device_launch_parameters.h>
#include <driver_types.h>
#include <cuda_runtime.h>
#include <stdio.h>

namespace NbodyCu {
void    Run(float timeStep, float4* X, float4* A, float4* V);
}
