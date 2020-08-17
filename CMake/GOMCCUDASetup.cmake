# Find CUDA is enabled, set it up

set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-DGOMC_CUDA)
set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -std=c++11")
include_directories(src/GPU)
#set(CUDA_VERBOSE_BUILD ON)
set(CUDA_SEPARABLE_COMPILATION ON)

set(GPU_NPT_flags "-DENSEMBLE=4 -DGOMC_CUDA")
set(GPU_NPT_name "GOMC_GPU_NPT")
set(GPU_GC_flags "-DENSEMBLE=3 -DGOMC_CUDA")
set(GPU_GC_name "GOMC_GPU_GCMC")
set(GPU_GE_flags "-DENSEMBLE=2 -DGOMC_CUDA")
set(GPU_GE_name "GOMC_GPU_GEMC")
set(GPU_NVT_flags "-DENSEMBLE=1 -DGOMC_CUDA")
set(GPU_NVT_name "GOMC_GPU_NVT")

#####################################
if(ENSEMBLE_GPU_NVT)
    cuda_add_executable(GPU_NVT ${cudaSources} ${cudaHeaders}
    ${sources} ${headers} ${libHeaders} ${libSources})
    set_target_properties(GPU_NVT PROPERTIES
        OUTPUT_NAME ${GPU_NVT_name}
        COMPILE_FLAGS "${GPU_NVT_flags}")
    if(WIN32)
        #needed for hostname
        target_link_libraries(GPU_NVT ws2_32)
    endif()
    if(MPI_FOUND)
	target_link_libraries(GPU_NVT ${MPI_LIBRARIES})
    endif()
endif()

if(ENSEMBLE_GPU_GEMC)
    cuda_add_executable(GPU_GEMC ${cudaSources} ${cudaHeaders} ${sources}
    ${headers} ${libHeaders} ${libSources})
    set_target_properties(GPU_GEMC PROPERTIES
        OUTPUT_NAME ${GPU_GE_name}
        COMPILE_FLAGS "${GPU_GE_flags}")
    if(WIN32)
        #needed for hostname
        target_link_libraries(GPU_GEMC ws2_32)
    endif()
    if(MPI_FOUND)
	target_link_libraries(GPU_GEMC ${MPI_LIBRARIES})
    endif()
endif()

if(ENSEMBLE_GPU_GCMC)
    cuda_add_executable(GPU_GCMC ${cudaSources} ${cudaHeaders} ${sources}
    ${headers} ${libHeaders} ${libSources})
    set_target_properties(GPU_GCMC PROPERTIES
        OUTPUT_NAME ${GPU_GC_name}
        COMPILE_FLAGS "${GPU_GC_flags}")
    if(WIN32)
        #needed for hostname
        target_link_libraries(GPU_GCMC ws2_32)
    endif()
    if(MPI_FOUND)
	target_link_libraries(GPU_GCMC ${MPI_LIBRARIES})
    endif()
endif()

if(ENSEMBLE_GPU_NPT)
    cuda_add_executable(GPU_NPT ${cudaSources} ${cudaHeaders} ${sources}
    ${headers} ${libHeaders} ${libSources})
    set_target_properties(GPU_NPT PROPERTIES
        OUTPUT_NAME ${GPU_NPT_name}
        COMPILE_FLAGS "${GPU_NPT_flags}")
    if(WIN32)
        #needed for hostname
        target_link_libraries(GPU_NPT ws2_32)
    endif()
    if(MPI_FOUND)
	target_link_libraries(GPU_NPT ${MPI_LIBRARIES})
    endif()
endif()
