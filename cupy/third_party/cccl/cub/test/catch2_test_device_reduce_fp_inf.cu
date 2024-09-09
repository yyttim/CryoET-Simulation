/******************************************************************************
 * Copyright (c) 2023, NVIDIA CORPORATION.  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the NVIDIA CORPORATION nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL NVIDIA CORPORATION BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 ******************************************************************************/

#include <cub/device/device_reduce.cuh>

#include <thrust/device_vector.h>

#include <cuda/std/limits>

// Has to go after all cub headers. Otherwise, this test won't catch unused
// variables in cub kernels.
#include "catch2/catch.hpp"
#include "catch2_test_cdp_helper.h"
#include "catch2_test_helper.h"

DECLARE_CDP_WRAPPER(cub::DeviceReduce::ArgMin, device_arg_min);
DECLARE_CDP_WRAPPER(cub::DeviceReduce::ArgMax, device_arg_max);

// %PARAM% TEST_CDP cdp 0:1

CUB_TEST("Device reduce arg{min,max} works with inf items", "[reduce][device]")
{
  using in_t     = float;
  using offset_t = int;
  using out_t    = cub::KeyValuePair<offset_t, in_t>;

  const int n     = 10;
  const float inf = ::cuda::std::numeric_limits<float>::infinity();

  thrust::device_vector<out_t> out(1);
  out_t *d_out = thrust::raw_pointer_cast(out.data());

  /**
   * ArgMin should return max value for empty input. This interferes with
   * input data containing infinity values. This test checks that ArgMin
   * works correctly with infinity values.
   */
  SECTION("InfInArgMin")
  {
    thrust::device_vector<in_t> in(n, inf);
    const in_t *d_in = thrust::raw_pointer_cast(in.data());

    device_arg_min(d_in, d_out, n);

    const out_t result = out[0];
    REQUIRE(result.key == 0);
    REQUIRE(result.value == inf);
  }

  /**
   * ArgMax should return lowest value for empty input. This interferes with
   * input data containing infinity values. This test checks that ArgMax
   * works correctly with infinity values.
   */
  SECTION("InfInArgMax")
  {
    thrust::device_vector<in_t> in(n, -inf);
    const in_t *d_in = thrust::raw_pointer_cast(in.data());

    device_arg_max(d_in, d_out, n);

    const out_t result = out[0];
    REQUIRE(result.key == 0);
    REQUIRE(result.value == -inf);
  }
}