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

#ifndef CHAPISUPERLU_H
#define CHAPISUPERLU_H

#include "chrono/core/ChPlatform.h"

// Chrono::Engine version
//
// This is an integer, as 0xaabbccdd where
// for example version 1.2.0 is 0x00010200

#define CH_VERSION_SUPERLU_MODULE 0x00000100

// When compiling this library, remember to define CH_API_COMPILE_SUPERLU
// (so that the symbols with 'ChApiSuperLU' in front of them will be
// marked as exported). Otherwise, just do not define it if you
// link the library to your code, and the symbols will be imported.

#if defined(CH_API_COMPILE_SUPERLU)
#define ChApiSuperLU ChApiEXPORT
#else
#define ChApiSuperLU ChApiIMPORT
#endif

/**
    @defgroup superlu_module SuperLU module
    @brief Module for the Intel MKL library direct solver

    Module provides access to the SuperLU library. This library is
    currently used in Chrono for its parallel direct solver (Pardiso).

    For additional information, see:
    - the [installation guide](@ref module_mkl_installation)
*/

#endif  // END of header
