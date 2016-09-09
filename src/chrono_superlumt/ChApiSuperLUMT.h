//
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2011-2012 Alessandro Tasora
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file at the top level of the distribution
// and at http://projectchrono.org/license-chrono.txt.
//

#ifndef CHAPISUPERLUMT_H
#define CHAPISUPERLUMT_H

#include "chrono/core/ChPlatform.h"

// Chrono::Engine version
//
// This is an integer, as 0xaabbccdd where
// for example version 1.2.0 is 0x00010200

#define CH_VERSION_SUPERLUT_MODULE 0x00000100

// When compiling this library, remember to define CH_API_COMPILE_SUPERLUMT
// (so that the symbols with 'ChApiSuperLUMT' in front of them will be
// marked as exported). Otherwise, just do not define it if you
// link the library to your code, and the symbols will be imported.

#if defined(CH_API_COMPILE_SUPERLUMT)
#define ChApiSuperLUMT ChApiEXPORT
#else
#define ChApiSuperLUMT ChApiIMPORT
#endif

/**
    @defgroup superlumt_module SuperLU_MT module
    @brief Module for the SuperLU_MT solver

    Module provides access to the SuperLU_MT library. This library is
    currently used in Chrono for its parallel direct solver (Pardiso).

    For additional information, see:
    - the [installation guide](@ref module_mkl_installation)
*/

#endif  // END of header
