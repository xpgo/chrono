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

#include "chrono/physics/ChSystem.h"
#include "chrono/utils/ChUtilsInputOutput.h"
#include "chrono_fea/ChElementShellTri_3.h"
#include "chrono_fea/ChLinkPointFrame.h"
#include "chrono_fea/ChMesh.h"
#include "chrono_irrlicht/ChIrrApp.h"

// basic file operations
#include <fstream>
#include <chrono_fea/ChVisualizationFEAmesh.h>

#define DUMP_LISTS true

using namespace chrono;
using namespace chrono::fea;
using namespace chrono::irrlicht;
using namespace irr;

int main(int argc, char* argv[]) {
    
    // Initialize
    auto my_mesh = std::make_shared<ChMesh>();

    ChSystem my_system;
    my_system.Set_G_acc(ChVector<>(0, 0, -9.8));

    // Create the Irrlicht visualization (open the Irrlicht device, bind a simple user interface, etc.)
    irrlicht::ChIrrApp application(&my_system, L"ANCF Shells", core::dimension2d<u32>(800, 600), false, true);

    // Easy shortcuts to add camera, lights, logo and sky in Irrlicht scene:
    application.AddTypicalLogo();
    application.AddTypicalSky();
    application.AddTypicalLights();
    application.AddTypicalCamera(core::vector3df(1,1,1),  // camera location
        core::vector3df(0.0f, 0.0f, 0.0f));  // "look at" location



    /****************************************
    * Beam creation
    ****************************************/
    // Input
    size_t cols_x = 3;
    size_t rows_y = 2;
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

            //if (col_sel == cols_x-1)
            //    node->SetPos(ChVector<>(loc_x_generator(col_sel), loc_y_generator(row_sel), -0.1));

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

    for (auto cont = 0; cont < my_mesh->GetNnodes(); cont++)
    {
        std::cout << "Node: " << my_mesh->GetNode(cont)->GetID() << std::endl;
        std::cout << "X0: " << std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(cont))->GetX0()(0);
        std::cout << "; X: " << std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(cont))->GetPos()(0) << std::endl;
        std::cout << "Y0: " << std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(cont))->GetX0()(1);
        std::cout << "; Y: " << std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(cont))->GetPos()(1) << std::endl;
        std::cout << "Z0: " << std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(cont))->GetX0()(2);
        std::cout << "; Z: " << std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(cont))->GetPos()(2) << std::endl;
        std::cout << std::endl;
    }

    // create Elements
    auto material = std::make_shared<ChMaterialShellTri_3>(210e9, 0.3, 7850);
    for (size_t col_sel = 0; col_sel < cols_x - 1; col_sel++)
    {
        for (size_t row_sel = 0; row_sel < rows_y - 1; row_sel++)
        {

            auto prova = std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(0));
            auto elementLEFTUP = std::make_shared<ChElementShellTri_3>();
            elementLEFTUP->SetNodes(std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode( row_sel     *cols_x + col_sel    )),
                                    std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode((row_sel + 1)*cols_x + col_sel + 1)),
                                    std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode((row_sel + 1)*cols_x + col_sel    )));

            elementLEFTUP->Set_thickness(0.1);
            elementLEFTUP->SetMaterial(material);
            my_mesh->AddElement(elementLEFTUP);

            auto elementBOTTOMRIGHT = std::make_shared<ChElementShellTri_3>();
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
            elem_list << std::dynamic_pointer_cast<ChElementGeneric>(my_mesh->GetElement(elem_sel))->GetID() << ", "
                << my_mesh->GetElement(elem_sel)->GetNodeN(0)->GetID() << ", "
                << my_mesh->GetElement(elem_sel)->GetNodeN(1)->GetID() << ", "
                << my_mesh->GetElement(elem_sel)->GetNodeN(2)->GetID() << std::endl;
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
        std::dynamic_pointer_cast<ChElementShellTri_3>(my_mesh->GetElement(elem_sel))->UpdateConnectivity(my_mesh);
    }

    // Add mesh to the system
    my_system.Add(my_mesh);
    my_system.SetupInitial();

    for (auto elem_sel = 0; elem_sel < my_mesh->GetNelements(); elem_sel++)
    {
        auto elem_temp = std::dynamic_pointer_cast<ChElementShellTri_3>(my_mesh->GetElement(elem_sel));
        std::cout << "Elem: " << elem_temp->GetID() << std::endl;
        std::cout << "Nodes: [";
        for (auto node_sel = 0; node_sel<6; ++node_sel)
        {
            if (elem_temp->all_nodes[node_sel])
                std::cout << elem_temp->all_nodes[node_sel]->GetID() << ", ";
            else
                std::cout << "np" << ", ";
        }
        std::cout << "\b\b]" << std::endl;
        std::cout << "Area0: " << elem_temp->GetArea0() << std::endl;
        std::cout << "EdgeLength0: " << "{" << elem_temp->edge_length0[0] << ", " << elem_temp->edge_length0[1] << ", " << elem_temp->edge_length0[2] << "}" << std::endl;
        std::cout << "EdgeVers0_0: " << "{" << elem_temp->edge_versors0[0](0) << ", " << elem_temp->edge_versors0[0](1) << ", " << elem_temp->edge_versors0[0](2) << "}" << std::endl;
        std::cout << "EdgeVers0_1: " << "{" << elem_temp->edge_versors0[1](0) << ", " << elem_temp->edge_versors0[1](1) << ", " << elem_temp->edge_versors0[1](2) << "}" << std::endl;
        std::cout << "EdgeVers0_2: " << "{" << elem_temp->edge_versors0[2](0) << ", " << elem_temp->edge_versors0[2](1) << ", " << elem_temp->edge_versors0[2](2) << "}" << std::endl;
        std::cout << "ElemNormal: " << "{" << elem_temp->edge_versors0[2](0) << ", " << elem_temp->edge_versors0[2](1) << ", " << elem_temp->edge_versors0[2](2) << "}" << std::endl;
        std::cout << std::endl;
    }


    

    // -------------------------------------
    // Irrlicht options
    // -------------------------------------

    auto mvisualizemesh = std::make_shared<ChVisualizationFEAmesh>(*(my_mesh.get()));
    mvisualizemesh->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NODE_SPEED_NORM);
    mvisualizemesh->SetColorscaleMinMax(0.0, 5.50);
    mvisualizemesh->SetShrinkElements(true, 0.85);
    mvisualizemesh->SetSmoothFaces(true);
    my_mesh->AddAsset(mvisualizemesh);

    auto mvisualizemeshref = std::make_shared<ChVisualizationFEAmesh>(*(my_mesh.get()));
    mvisualizemeshref->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_SURFACE);
    mvisualizemeshref->SetWireframe(true);
    mvisualizemeshref->SetDrawInUndeformedReference(true);
    my_mesh->AddAsset(mvisualizemeshref);

    auto mvisualizemeshC = std::make_shared<ChVisualizationFEAmesh>(*(my_mesh.get()));
    mvisualizemeshC->SetFEMglyphType(ChVisualizationFEAmesh::E_GLYPH_NODE_DOT_POS);
    mvisualizemeshC->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NONE);
    mvisualizemeshC->SetSymbolsThickness(0.004);
    my_mesh->AddAsset(mvisualizemeshC);

    auto mvisualizemeshD = std::make_shared<ChVisualizationFEAmesh>(*(my_mesh.get()));
    // mvisualizemeshD->SetFEMglyphType(ChVisualizationFEAmesh::E_GLYPH_NODE_VECT_SPEED);
    mvisualizemeshD->SetFEMglyphType(ChVisualizationFEAmesh::E_GLYPH_ELEM_TENS_STRAIN);
    mvisualizemeshD->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NONE);
    mvisualizemeshD->SetSymbolsScale(1);
    mvisualizemeshD->SetColorscaleMinMax(-0.5, 5);
    mvisualizemeshD->SetZbufferHide(false);
    my_mesh->AddAsset(mvisualizemeshD);

    application.AssetBindAll();
    application.AssetUpdateAll();
    
    // Set up integrator
    my_system.SetIntegrationType(ChSystem::INT_HHT);
    auto mystepper = std::static_pointer_cast<ChTimestepperHHT>(my_system.GetTimestepper());
    mystepper->SetAlpha(-0.2);
    mystepper->SetMaxiters(100);
    mystepper->SetAbsTolerances(1e-5);
    mystepper->SetMode(ChTimestepperHHT::POSITION);
    mystepper->SetScaling(true);
    mystepper->SetVerbose(true);


    ChMatrixDynamic<double> H;
    ChMatrixDynamic<double> Fi;
    my_mesh->GetElement(0)->ComputeInternalForces(Fi);
    my_mesh->GetElement(0)->ComputeKRMmatricesGlobal(H, 1, 0, 0);

    GetLog() << H << "\n";
    GetLog() << Fi;

    //int num_steps = 2;
    //double step_size = 0.01;
    //for (int istep = 0; istep < num_steps; istep++) {
    //    my_system.DoStepDynamics(step_size);
    //}

    application.SetTimestep(0.001);
    application.BeginScene();
    application.DrawAll();
    application.EndScene();

    while (application.GetDevice()->run()) {
        application.BeginScene();
        application.DrawAll();
        application.DoStep();
        application.EndScene();
    }

    
    return 0;
}



