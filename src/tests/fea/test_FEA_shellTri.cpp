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



    /****************************************
    * Beam creation
    ****************************************/
    // Input
    int cols_x = 3;
    int rows_y = 2;
    double length = 2;
    double width = 1;

    // Start creation
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

    // create nodes
    for (size_t row_sel = 0; row_sel < rows_y; row_sel++)
    {
        for (size_t col_sel = 0; col_sel < cols_x; col_sel++)
        {
            auto node = std::make_shared<ChNodeFEAxyz>(ChVector<>(loc_x_generator(col_sel), loc_y_generator(row_sel), 0));

            if (col_sel == 0) // fix the base nodes
                node->SetFixed(true);

            my_mesh->AddNode(node);
        }
    }

    
    if (DUMP_LISTS)
    {
        std::ofstream nodes_list;
        nodes_list.open("nodes_list.txt");
        nodes_list << "# Nodes list[NODEID, X, Y, Z]" << std::endl;
        for (size_t node_sel = 0; node_sel < my_mesh->GetNnodes(); ++node_sel)
        {
            nodes_list << my_mesh->GetNode(node_sel).get() << ", "
                << std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(node_sel))->GetPos().x << ", "
                << std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(node_sel))->GetPos().y << ", "
                << std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(node_sel))->GetPos().z << std::endl;
        }
        nodes_list.close();
    }



    // create Elements
    auto material = std::make_shared<ChMaterialShellTri>(210e9, 0.3, 7850);
    for (size_t col_sel = 0; col_sel < cols_x - 1; col_sel++)
    {
        for (size_t row_sel = 0; row_sel < rows_y - 1; row_sel++)
        {

            auto prova = std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(0));
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
                                         std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode( row_sel   *cols_x + col_sel+1 )));

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
                << my_mesh->GetElement(elem_sel)->GetNodeN(0).get() << ", "
                << my_mesh->GetElement(elem_sel)->GetNodeN(1).get() << ", "
                << my_mesh->GetElement(elem_sel)->GetNodeN(2).get() << std::endl;
        }

        elem_list.close();
    }

    // Switch off mesh class gravity
    my_mesh->SetAutomaticGravity(false);

    // Let element find its neighbour
    // TODO: finding the neighbours is mandatory; it could be done inside the SetupInitial()...
    // TODO: .. but the SetupInitial() member function overrides a base function with no arguments ...
    // TODO: ... instead I need to pass the ChMesh that contains a list of other elements within which search the neighbours
    for (size_t elem_sel = 0; elem_sel < my_mesh->GetNelements(); elem_sel++)
    {
        std::dynamic_pointer_cast<ChElementShellTri>(my_mesh->GetElement(elem_sel))->UpdateConnectivity(my_mesh);
    }


    // Add mesh to the system
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
        my_system.DoStepDynamics(step_size);
    }

    
    return 0;
}
