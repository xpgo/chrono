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
#include "chrono/solver/ChSolverMINRES.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChSystem.h"
#include "chrono/utils/ChUtilsInputOutput.h"
#include "chrono_fea/ChElementShellTri.h"
#include "chrono_fea/ChLinkDirFrame.h"
#include "chrono_fea/ChLinkPointFrame.h"
#include "chrono_fea/ChMesh.h"

// basic file operations
#include <iostream>
#include <fstream>
#include <chrono_fea/ChElementShellANCF.h>

#ifdef CHRONO_MKL
#include "chrono_mkl/ChSolverMKL.h"
////#define USE_MKL
#else
#undef USE_MKL
#endif

#ifdef CHRONO_OPENMP_ENABLED
#include <omp.h>
#endif

#define DUMP_LISTS true

using namespace chrono;
using namespace chrono::fea;


int main(int argc, char* argv[]) {
    
    // Initialize
    auto my_mesh = std::make_shared<ChMesh>();
    std::vector<ChElementShellTri> TriShell;

    // Input
    int cols_x = 3;
    int rows_y = 2;
    double length = 1;
    double width = 0.1;

    // Beam creation
    double x_step = length / (cols_x - 1);
    double y_step = width / (rows_y - 1);

    chrono::ChVectorDynamic<double> loc_x_generator(rows_y);
    chrono::ChVectorDynamic<double> loc_y_generator(cols_x);
    loc_x_generator.Resize(cols_x);
    loc_y_generator.Resize(rows_y);

    // x coordinates
    for (size_t col_sel = 0; col_sel < cols_x; col_sel++)
    {
        loc_x_generator(col_sel) = col_sel*x_step;
    }

    // y coordinates
    for (size_t row_sel = 0; row_sel < rows_y; row_sel++)
    {
        loc_y_generator(row_sel) = row_sel*y_step;
    }

    std::ofstream nodes_list;
    if (DUMP_LISTS)
    {
        nodes_list.open("nodes_list.txt");
        nodes_list << "# Nodes list[NODEID, X, Y, Z]" << std::endl;
    }
    
    

    // create nodes
    for (size_t col_sel = 0; col_sel < cols_x; col_sel++)
    {
        for (size_t row_sel = 0; row_sel < rows_y; row_sel++)
        {
            auto node = std::make_shared<chrono::fea::ChNodeFEAxyz>(chrono::ChVector<>(loc_x_generator(col_sel), loc_y_generator(row_sel), 0));
            my_mesh->AddNode(node);

            if (DUMP_LISTS)
            {
                nodes_list << node->NodeGetOffset_w() << ", "
                    << node->GetPos()(0) << ", "
                    << node->GetPos()(1) << ", "
                    << node->GetPos()(2) << std::endl;
            }
            

            if (row_sel == 0) // fix the base nodes
                node->SetFixed(true);
        }
    }

    

    // create Elements
    auto material = std::make_shared<ChMaterialShellTri>(210e9, 0.3, 7850);
    for (size_t col_sel = 0; col_sel < cols_x - 1; col_sel++)
    {
        for (size_t row_sel = 0; row_sel < rows_y - 1; row_sel++)
        {
            auto elementLEFTUP = std::make_shared<ChElementShellTri>();
            elementLEFTUP->SetNodes(std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode( row_sel     *cols_x + col_sel    )),
                                    std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode((row_sel + 1)*cols_x + col_sel + 1)),
                                    std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode((row_sel + 1)*cols_x + col_sel    )));

            elementLEFTUP->Set_thickness(0.1);
            elementLEFTUP->SetMaterial(material);
            my_mesh->AddElement(elementLEFTUP);

            auto elementBOTTOMRIGHT = std::make_shared<ChElementShellTri>();
            elementBOTTOMRIGHT->SetNodes(std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode( row_sel   *cols_x + col_sel   )),
                                         std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode((row_sel+1)*cols_x + col_sel+1 )),
                                         std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(row_sel    *cols_x + col_sel+1 )));

            elementBOTTOMRIGHT->Set_thickness(0.1);
            elementBOTTOMRIGHT->SetMaterial(material);
            my_mesh->AddElement(elementBOTTOMRIGHT);
        }
    }


    if (DUMP_LISTS)
    {
        std::ofstream elem_list;

        elem_list.open("elem_list.txt");
        elem_list << "# Element list [ELEMENTID, NODE1, NODE2, NODE3]" << std::endl;
        for (size_t elem_sel = 0; elem_sel < my_mesh->GetNelements(); ++elem_sel)
        {
            elem_list << elem_sel << ", "
                << my_mesh->GetNode(0)->NodeGetOffset_w() << ", "
                << my_mesh->GetNode(1)->NodeGetOffset_w() << ", "
                << my_mesh->GetNode(2)->NodeGetOffset_w() << std::endl;
        }

        elem_list.close();
    }

    // Switch off mesh class gravity (ANCF shell elements have a custom implementation)
    my_mesh->SetAutomaticGravity(false);

    for (size_t elem_sel = 0; elem_sel < my_mesh->GetNelements(); elem_sel++)
    {
        std::dynamic_pointer_cast<ChElementShellTri>(my_mesh->GetElement(elem_sel))->UpdateConnectivity(my_mesh);
    }


    ChSystem my_system;
    my_system.Add(my_mesh);
    my_system.SetupInitial();

    
    // Set up integrator
    my_system.SetIntegrationType(ChSystem::INT_HHT);
    auto mystepper = std::static_pointer_cast<ChTimestepperHHT>(my_system.GetTimestepper());
    mystepper->SetAlpha(-0.2);
    mystepper->SetMaxiters(100);
    mystepper->SetAbsTolerances(1e-5);
    mystepper->SetMode(ChTimestepperHHT::POSITION);
    mystepper->SetScaling(true);
    //// mystepper->SetVerbose(true);

    int num_steps = 2;
    double step_size = 0.01;
    for (int istep = 0; istep < num_steps; istep++) {
        ////GetLog() << " step number: " << istep << "  time: " << my_system.GetChTime() << "\n";
        my_system.DoStepKinematics(step_size);
    }

    
    return 0;
}
