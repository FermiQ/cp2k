/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2025 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#ifndef DBM_MULTIPLY_GPU_KERNEL_H
#define DBM_MULTIPLY_GPU_KERNEL_H

#include "../offload/offload_runtime.h"
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_DBM)

#include "dbm_internal.h"

#ifdef __cplusplus
extern "C" {
#endif

/*******************************************************************************
 * \brief Internal routine for launching the GPU kernel.
 *        All arguments are assumed to be device pointers.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_multiply_gpu_launch_kernel(offloadStream_t stream, double alpha,
                                    int ntasks, const dbm_task_t *tasks_host,
                                    const dbm_task_t *tasks,
                                    const double *pack_a_data,
                                    const double *pack_b_data,
                                    double *shard_c_data);

#ifdef __cplusplus
}
#endif

#endif // defined(__OFFLOAD) && !defined(__NO_OFFLOAD_DBM)
#endif

// EOF
