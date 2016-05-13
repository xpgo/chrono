//
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2013 Project Chrono
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file at the top level of the distribution
// and at http://projectchrono.org/license-chrono.txt.
//

#include "chrono/ChConfig.h"
#include "chrono/core/ChFileutils.h"
#include "chrono/core/ChTimer.h"
#include "chrono/lcp/ChLcpIterativeMINRES.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChSystem.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#include "chrono_fea/ChElementShellTri.h"
#include "chrono_fea/ChLinkDirFrame.h"
#include "chrono_fea/ChLinkPointFrame.h"
#include "chrono_fea/ChMesh.h"

#ifdef CHRONO_MKL
#include "chrono_mkl/ChLcpMklSolver.h"
////#define USE_MKL
#else
#undef USE_MKL
#endif

#ifdef CHRONO_OPENMP_ENABLED
#include <omp.h>
#endif

void beam_custom_creator( std::shared_ptr<ChMesh> my_mesh, size_t cols_x, size_t rows_y, double length, double width ) {

    double x_step = length / (cols_x - 1);
    double y_step = width / (rows_y - 1);

    ChVectorDynamic<double> loc_x_generator(rows_y);
    ChVectorDynamic<double> loc_y_generator(cols_x);

    for (size_t col_sel = 0; col_sel < cols_x; col_sel++)
    {
        loc_x_generator(col_sel) = col_sel*x_step;
    }

    for (size_t row_sel = 0; row_sel < rows_y; row_sel++)
    {
        loc_y_generator(row_sel) = row_sel*y_step;
    }

    for (size_t col_sel = 0; col_sel < cols_x; col_sel++)
    {
        for (size_t row_sel = 0; row_sel < rows_y; row_sel++)
        {
            auto node = std::make_shared<ChNodeFEAxyz>(ChVector<>(loc_x_generator(col_sel), loc_y_generator(row_sel), 0));
            my_mesh->AddNode(node);

            if (row_sel == 0) // fix the base nodes
                node->SetFixed(true);
        }
    }

    
}

using namespace chrono;
using namespace chrono::fea;


int main(int argc, char* argv[]) {
    
    auto my_mesh = std::make_shared<ChMesh>();
    
    return 0;
}
