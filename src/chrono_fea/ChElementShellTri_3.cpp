#include "chrono_fea/ChElementShellTri_3.h"
#include <bitset>

namespace chrono {
namespace fea {

int ChElementShellTri_3::getLocalNodeNumber(int element, int neighbour_numbering)
{
    return neigh_to_local_numbering[element][neighbour_numbering];
}

    int ChElementShellTri_3::countNeighbours() const
    {
        int count = 0;
        for (auto node_sel = 3; node_sel<6; ++node_sel)
        {
            if (all_nodes[node_sel])
            {
                count++;
            }
        }
        return count;
    }

    void ChElementShellTri_3::updateBC()
    {
        if (all_nodes[0]->GetFixed() && all_nodes[1]->GetFixed())
            edge_bc[2] = boundary_conditions::CLAMPED;

        if (all_nodes[1]->GetFixed() && all_nodes[2]->GetFixed())
            edge_bc[0] = boundary_conditions::CLAMPED;

        if (all_nodes[2]->GetFixed() && all_nodes[0]->GetFixed())
            edge_bc[1] = boundary_conditions::CLAMPED;
    }

    void ChElementShellTri_3::updateElementMass(int node_sel)
    {
        mass = element_area0 * thickness * m_material->Get_rho();
        for (auto diag_sel = 0; diag_sel < 3; diag_sel++)
        {
            all_nodes[diag_sel]->SetMass(mass / 3);
        }
    }

    void ChElementShellTri_3::getNodeM(ChMatrix<double>& node_mass, int node_sel, double factor)
    {
        if (node_mass.GetColumns() != 3 || node_mass.GetRows() != 3)
            node_mass.Resize(3, 3);

        updateElementMass(node_sel);
        node_mass.FillDiag(all_nodes[node_sel]->GetMass() * factor);
    }

    void ChElementShellTri_3::getNodeK(ChMatrix<double>& node_damp, int node_sel, double factor)
    {
        if (node_damp.GetColumns() != 3 || node_damp.GetRows() != 3)
            node_damp.Resize(3, 3);

        node_damp.FillDiag(1e-3);
    }

    void ChElementShellTri_3::getElementMR(ChMatrix<>& H, double Mfactor, double Kfactor)
    {
        H.Resize(GetNdofs(), GetNdofs());

        //// Fill the H matrix with damping and mass; K and M supposed diagonal-block
        //ChMatrix33<double> node_mat;
        //for (auto main_node_sel = 0; main_node_sel<3; main_node_sel++)
        //{
        //    getNodeK(node_mat, main_node_sel, Kfactor);
        //    H.PasteSumMatrix(&node_mat, main_node_sel * 3, main_node_sel * 3);
        //    getNodeM(node_mat, main_node_sel, Mfactor);
        //    H.PasteSumMatrix(&node_mat, main_node_sel * 3, main_node_sel * 3);
        //}

        //int offset_diag = 0;
        //for (auto diag_sel = 0; diag_sel<3; ++diag_sel)
        //{
        //    if (neighbouring_elements[diag_sel].get() == nullptr)
        //    {
        //        offset_diag++;
        //        continue;
        //    }

        //    neighbouring_elements[diag_sel]->getNodeK(node_mat, neighbour_node_not_shared[diag_sel], Kfactor);
        //    H.PasteSumMatrix(&node_mat, (diag_sel + 3 - offset_diag) * 3, (diag_sel + 3 - offset_diag) * 3);
        //    neighbouring_elements[diag_sel]->getNodeM(node_mat, neighbour_node_not_shared[diag_sel], Mfactor);
        //    H.PasteSumMatrix(&node_mat, (diag_sel + 3 - offset_diag) * 3, (diag_sel + 3 - offset_diag) * 3);
        //}
    }

    void ChElementShellTri_3::getElementData(const ChVector<>& P1_in, const ChVector<>& P2_in, const ChVector<>& P3_in,
                                             std::array<double, 3>* edge_length,
                                             ChVector<>* z_versor,
                                             double* area,
                                             ChMatrix33<>* rotGL,
                                             ChMatrix<>* gradient_shape_function,
                                             ChMatrix<>* gradient_local,
                                             ChMatrix<>* gradient_side_n,
                                             ChMatrix<>* edge_normal_vers,
                                             ChMatrix<>* c_proj
                                             )
    {
        double area_temp;
        if (!area)
            area = &area_temp;

        // edge vectors
        ChVector<double> e_vect1;
        ChVector<double> e_vect2;
        ChVector<double> e_vect3;

        e_vect1.Sub(P3_in, P2_in);
        e_vect2.Sub(P1_in, P3_in);
        e_vect3.Sub(P2_in, P1_in);

        if (edge_length || c_proj)
        {
            (*edge_length)[0] = e_vect1.Length();
            (*edge_length)[1] = e_vect2.Length();
            (*edge_length)[2] = e_vect3.Length();
        }
        
        // local Cartesian axis
        ChVector<double> x_versor;
        ChVector<double> y_versor;
        
        z_versor->Cross(e_vect3, -e_vect2);
        auto lambda = z_versor->Length();
        *area = lambda / 2.0f; //TODO: strange fact...
        z_versor->Scale(1.0 / lambda);
        x_versor = e_vect3;
        x_versor.Normalize();
        y_versor.Cross(*z_versor, x_versor);

        if (rotGL || gradient_shape_function || gradient_local)
        {
            // rotation matrix from Local coordinates to Global coordinates
            rotGL->PasteVector(x_versor, 0, 0);
            rotGL->PasteVector(y_versor, 0, 1);
            rotGL->PasteVector(*z_versor, 0, 2);
        }
        

        if (gradient_shape_function || gradient_local)
        {
            ChMatrix33<double> rotLG;
            rotGL->FastInvert(&rotLG); //TODO: can be done in-place on rotGL?

            ChMatrix33<double> nodes_local;
            //nodes_local.PasteVector(rotLG.Matr_x_Vect(P1_in - P1_in) , 0, 0); 
            nodes_local.PasteVector(rotLG.Matr_x_Vect(P2_in - P1_in), 0, 1);
            nodes_local.PasteVector(rotLG.Matr_x_Vect(P3_in - P1_in), 0, 2);

            // normal-to-edge vectors
            ChMatrixNM<double, 2, 3> edge_normal_vect;
            edge_normal_vect(0, 0) = +nodes_local(1, 2) - nodes_local(1, 1);
            edge_normal_vect(0, 1) = +nodes_local(1, 0) - nodes_local(1, 2);
            edge_normal_vect(0, 2) = +nodes_local(1, 1) - nodes_local(1, 0);
            edge_normal_vect(1, 0) = -nodes_local(0, 2) + nodes_local(0, 1);
            edge_normal_vect(1, 1) = -nodes_local(0, 0) + nodes_local(0, 2);
            edge_normal_vect(1, 2) = -nodes_local(0, 1) + nodes_local(0, 0);

            // gradient shape function
            gradient_shape_function->CopyFromMatrixT(edge_normal_vect);
            *gradient_shape_function *= -0.5f / *area;
        }
        

        if (gradient_local || gradient_side_n)
        {
            // gradients (local frame)
            ChMatrix33<double> nodes_glob;
            nodes_glob.PasteVector(P1_in, 0, 0);
            nodes_glob.PasteVector(P2_in, 0, 1);
            nodes_glob.PasteVector(P3_in, 0, 2);
            gradient_local->MatrMultiply(nodes_glob, *gradient_shape_function);
        }
        

        if (edge_normal_vers || gradient_side_n)
        {
            double norm_temp;
            norm_temp = sqrt(pow((*edge_normal_vers)(0, 0), 2) + pow((*edge_normal_vers)(1, 0), 2));
            (*edge_normal_vers)(0, 0) /= norm_temp;
            (*edge_normal_vers)(1, 0) /= norm_temp;
            norm_temp = sqrt(pow((*edge_normal_vers)(0, 1), 2) + pow((*edge_normal_vers)(1, 1), 2));
            (*edge_normal_vers)(0, 1) /= norm_temp;
            (*edge_normal_vers)(1, 1) /= norm_temp;
            norm_temp = sqrt(pow((*edge_normal_vers)(0, 2), 2) + pow((*edge_normal_vers)(1, 2), 2));
            (*edge_normal_vers)(0, 2) /= norm_temp;
            (*edge_normal_vers)(1, 2) /= norm_temp;
        }
        

        if (gradient_side_n)
        {
            // gradients (side frame)
            ChMatrixNM<double, 2, 1> edge_normal_vers_temp;
            edge_normal_vers_temp(0, 0) = -(*edge_normal_vers)(0, 2);
            edge_normal_vers_temp(1, 0) = -(*edge_normal_vers)(1, 2);
            gradient_side_n->MatrMultiply(*gradient_local, edge_normal_vers_temp);
        }

        if (c_proj)
        {
            // C projections
            ChMatrix33<double> c_temp1;
            c_temp1.PasteVector(e_vect1, 0, 0);
            c_temp1.PasteVector(e_vect2, 0, 1);
            c_temp1.PasteVector(e_vect3, 0, 2);

            ChMatrix33<double> c_temp2;
            c_temp2.PasteVector(e_vect1 / (*edge_length)[0], 0, 0);
            c_temp2.PasteVector(e_vect2 / (*edge_length)[1], 0, 1);
            c_temp2.PasteVector(e_vect3 / (*edge_length)[2], 0, 2);

            c_proj->MatrTMultiply(c_temp1, c_temp2);
            *c_proj *= 0.5f / *area;
        }
        
        

    }

    void ChElementShellTri_3::Update()
    {
        ChMatrixNM<double,18,1> lambda_gamma;
        ChMatrixNM<double, 3, 18> Bb;
        ChMatrixNM<double, 3, 18> Bm;

        ChVector<double> z_versor;
        std::array<double, 3> edge_length = { -1,-1,-1 };
        ChMatrix33<double> rotGL;
        ChMatrixNM<double, 2, 3> edge_normal_vers;
        ChMatrixNM<double, 3, 2> gradient_shape_function;
        ChMatrixNM<double, 3, 2> gradient_local;
        ChMatrixNM<double, 3, 1> gradient_side_n;
        ChMatrix33<double> c_proj;
        double element_area = -1;

        // variables that hold values for neighbours
        ChVector<double> z_versor_adj;
        std::array<double, 3> edge_length_adj = { -1,-1,-1 };
        ChMatrix33<double> rotGL_adj;
        ChMatrixNM<double, 2, 3> edge_normal_vers_adj;
        ChMatrixNM<double, 3, 2> gradient_shape_function_adj;
        ChMatrixNM<double, 3, 2> gradient_local_adj;
        ChMatrixNM<double, 3, 1> gradient_side_n_adj;
        ChMatrix33<double> c_proj_adj;
        double element_area_adj = -1;

        getElementData(all_nodes[0]->GetPos(), all_nodes[1]->GetPos(), all_nodes[2]->GetPos(),
                       &edge_length,
                       &z_versor,
                       &element_area,
                       &rotGL,
                       &gradient_shape_function,
                       &gradient_local,
                       &gradient_side_n,
                       &edge_normal_vers,
                       &c_proj );

        update_counter++;
        

        for (auto edge_sel = 0; edge_sel<3; edge_sel++)
        {
            //TODO: should be done only once per Update()

            getElementData(all_nodes[getLocalNodeNumber(edge_sel + 1, 0)]->GetPos(),
                           all_nodes[getLocalNodeNumber(edge_sel + 1, 1)]->GetPos(),
                           all_nodes[getLocalNodeNumber(edge_sel + 1, 2)]->GetPos(),
                           &edge_length_adj,
                           &z_versor_adj,
                           &element_area_adj,
                           &rotGL_adj,
                           &gradient_shape_function_adj,
                           &gradient_local_adj,
                           &gradient_side_n_adj,
                           &edge_normal_vers_adj,
                           &c_proj_adj);


            // compute lamba*gamma
            lambda_gamma.FillElem(0);
            lambda_gamma.PasteVector(z_versor*c_proj0(0,edge_sel), 0, 0);
            lambda_gamma.PasteVector(z_versor*c_proj0(1,edge_sel), 3, 0);
            lambda_gamma.PasteVector(z_versor*c_proj0(2,edge_sel), 6, 0);

            if (neighbouring_elements[edge_sel]) //TODO: the numbering of neighbours is not predictable: next line must be changed accordingly
            {
                for (auto node_sel = 0; node_sel < 3; node_sel++)
                {
                    lambda_gamma.PasteSumVector(z_versor_adj*neighbouring_elements[edge_sel]->c_proj0(node_sel, 2),
                                                getLocalNodeNumber(edge_sel + 1,node_sel)*3, 0);

                }
            }

            

            // relative stiffness
            // TODO: relate rM to Young mod and thickness
            double riM = 0.5f;

            // curvature matrix
            
            ChMatrixNM<double, 3, 18> Bb_temp;
            Bb_temp.MatrMultiplyT(L_block0[edge_sel], lambda_gamma);
            Bb_temp *= riM / edge_length0[edge_sel];
            Bb.PasteSumMatrix(&Bb_temp,0,0);

            ChMatrixNM<double, 3, 3> Bm_temp;
            Bm_temp.PasteVector(gradient_local.ClipVector(0, 0)*gradient_shape_function0(edge_sel, 0), 0, 0);
            Bm_temp.PasteVector(gradient_local.ClipVector(0, 1)*gradient_shape_function0(edge_sel, 1), 0, 1);
            Bm_temp.PasteVector(gradient_local.ClipVector(0, 1)*gradient_shape_function0(edge_sel, 0) +
                                gradient_local.ClipVector(0, 0)*gradient_shape_function0(edge_sel, 1), 0, 2);
            Bm.PasteTranspMatrix(&Bm_temp, 0, edge_sel * 3);
        }

        
        Bb *= 4.0f * element_area0;
        
        // stiffness matrix
        // membrane
        ChMatrixNM<double, 3, 18> temp;
        temp.MatrMultiply(m_material->GetConsitutiveMatrix(), Bm);
        temp *= thickness0;
        stiffness_matrix.MatrTMultiply(Bm, temp);
        // bending
        temp.MatrMultiply(m_material->GetConsitutiveMatrix(), Bb);
        temp *= pow(thickness0,3)/12.0f;
        ChMatrixNM<double, 18, 18> stiffness_matrix_temp;
        stiffness_matrix_temp.MatrTMultiply(Bb, temp);
        stiffness_matrix += stiffness_matrix_temp;
        stiffness_matrix *= element_area0;

    }



    void ChElementShellTri_3::UpdateConnectivity(std::shared_ptr<ChMesh> mesh)
    {
        for (unsigned int elem_sel = 0; elem_sel < mesh->GetNelements(); ++elem_sel) // pick an element in the list that is candidate to be a neighbour
        {
            std::bitset<3> main_node_shared;
            std::shared_ptr<ChNodeFEAbase> not_shared_node_ptr;
            for (auto neigh_node_sel = 0; neigh_node_sel < 3; ++neigh_node_sel) // pick a node of the candidate element
            {
                int main_node_sel;
                for (main_node_sel = 0; main_node_sel < 3; ++main_node_sel) // pick a node of the current element
                {
                    if (mesh->GetElement(elem_sel)->GetNodeN(neigh_node_sel) == all_nodes[main_node_sel]) // if the two nodes picked are different...
                    {
                        main_node_shared[main_node_sel] = true;
                        break;
                    }
                }

                if (main_node_sel == 3)
                {
                    if (!not_shared_node_ptr)
                    {
                        not_shared_node_ptr = mesh->GetElement(elem_sel)->GetNodeN(neigh_node_sel);
                    }
                    else
                        break;
                }
            }

            if (not_shared_node_ptr && main_node_shared.count() == 2) // exclude the case that the candidate and current elements coincide;
            {
                // on which edge the just-found neighbour is attached?
                int not_shared_node_sel; // the not_shared_node_sel will point to the node in all_nodes that is not shared
                // --> it will be also the number of the edge on which it is attached
                for (not_shared_node_sel = 0; not_shared_node_sel < 3; ++not_shared_node_sel)
                {
                    if (!main_node_shared.test(not_shared_node_sel))
                    {
                        break;
                    }
                }

                all_nodes[3 + not_shared_node_sel] = std::dynamic_pointer_cast<ChNodeFEAxyz>(not_shared_node_ptr);
                neighbouring_elements[not_shared_node_sel] = std::dynamic_pointer_cast<ChElementShellTri_3>(mesh->GetElement(elem_sel));
            }
        }
    }

    void ChElementShellTri_3::SetBC(int BC_on_edge1, int BC_on_edge2, int BC_on_edge3)
    {
        edge_bc[0] = static_cast<boundary_conditions>(BC_on_edge1);
        edge_bc[1] = static_cast<boundary_conditions>(BC_on_edge2);
        edge_bc[2] = static_cast<boundary_conditions>(BC_on_edge3);
    }

    void ChElementShellTri_3::ComputeKRMmatricesGlobal(ChMatrix<>& H, double Kfactor, double Rfactor, double Mfactor)
    {
        // WARNING: the stiffness matrix is supposed to be already updated since ComputeInternalForces()
        // should have been called previously

        H.Resize(GetNdofs(), GetNdofs());

        // the N-W corner is always present
        H.PasteClippedMatrix(&stiffness_matrix, 0, 0, 9, 9, 0, 0);

        // TODO: check if the stiffness matrix is symmetric? I think it isn't...
        int offset_row = 0;
        int offset_col;
        for (auto row_sel = 0; row_sel < 3; row_sel++)
        {
            if (neighbouring_elements[row_sel].get() == nullptr)
            {
                offset_row++;
                continue;
            }

            offset_col = 0;
            for (auto col_sel = 0; col_sel < 3; col_sel++)
            {
                if (neighbouring_elements[row_sel].get() == nullptr)
                {
                    offset_col++;
                    continue;
                }

                H.PasteClippedMatrix(&stiffness_matrix, row_sel * 3, col_sel * 3, 3, 3, (row_sel - offset_row) * 3, (col_sel - offset_col) * 3);
            }
        }

        H.MatrScale(Kfactor);

        getElementMR(H, Mfactor, Kfactor);
    }

    void ChElementShellTri_3::ComputeInternalForces(ChMatrixDynamic<>& Fi)
    {

        // TODO: it shouldn't be necessary since K it is updated during Update() call
        Update();

        ChMatrixNM<double, 18, 1> disp_glob;
        // paste main_nodes position vectors
        for (auto node_sel = 0; node_sel < 6; node_sel++)
        {
            if (!all_nodes[node_sel])
            {
                disp_glob(node_sel * 3) = 0;
                disp_glob(node_sel * 3 + 1) = 0;
                disp_glob(node_sel * 3 + 2) = 0;
            }
            else
                disp_glob.PasteVector(all_nodes[node_sel]->GetPos() - all_nodes[node_sel]->GetX0(), node_sel * 3, 0);
        }

        
        if (GetNdofs() == 18)
        {
            Fi.Resize(18, 1);
            Fi.MatrMultiply(stiffness_matrix, disp_glob);
        }
        else
        {
            Fi.Resize(GetNdofs(), 1);

            ChMatrixNM<double, 18, 1> Fi_temp;
            Fi_temp.MatrMultiply(stiffness_matrix, disp_glob);

            // paste forces of the three main nodes
            Fi.PasteClippedMatrix(&Fi_temp, 0, 0, 9, 1, 0, 0);

            // paste forces of the present neigh nodes
            int offset_row = 0;
            for (auto neigh_elem_sel = 0; neigh_elem_sel < 3; neigh_elem_sel++)
            {
                if (neighbouring_elements[neigh_elem_sel])
                {
                    Fi.PasteClippedMatrix(&Fi_temp, (neigh_elem_sel + 3) * 3, 0, 3, 1, (neigh_elem_sel + 3 - offset_row) * 3, 0);
                }
                else
                    offset_row++;
            }
        }
    }

    void ChElementShellTri_3::SetupInitial(ChSystem* system)
    {
        m_material->UpdateConsitutiveMatrices();
        updateBC();

        getElementData(all_nodes[0]->GetX0(), all_nodes[1]->GetX0(), all_nodes[2]->GetX0(),
                       &edge_length0,
                       &z_versor0,
                       &element_area0,
                       &rotGL0,
                       &gradient_shape_function0,
                       &gradient_local0,
                       &gradient_side_n0,
                       &edge_normal_vers0,
                       &c_proj0);

        // compute L block
        for (auto edge_sel = 0; edge_sel<3; ++edge_sel)
        {
            L_block0[edge_sel](0) = pow(gradient_shape_function0(edge_sel, 0), 2);
            L_block0[edge_sel](1) = pow(gradient_shape_function0(edge_sel, 1), 2);
            L_block0[edge_sel](2) = 2 * gradient_shape_function0(edge_sel, 0)*gradient_shape_function0(edge_sel, 1);
        }

        // Inform the system about which nodes take part to the computation
        //TODO: do I have to put also neighbouring nodes?
        std::vector<ChVariables*> vars;
        for (auto node_sel : all_nodes)
        {
            if (node_sel)
                vars.push_back(&node_sel->Variables());
        }

        Kmatr.SetVariables(vars);
    }

    std::shared_ptr<ChNodeFEAbase> ChElementShellTri_3::GetNodeN(int n)
    {
        int offset = 0;
        for (auto node_sel = 0; node_sel <= n + offset; node_sel++)
        {
            if (!all_nodes[node_sel])
                offset++;
        }

        return all_nodes[n + offset];
    }

    void ChElementShellTri_3::GetStateBlock(ChMatrixDynamic<>& mD)
    {
        mD.Resize(GetNdofs(), 1);
        int offset_row = 0;
        for (auto node_sel : all_nodes)
        {
            if (node_sel)
            {
                mD.PasteVector(node_sel->GetPos(), offset_row, 0);
                offset_row += 3;
            }
        }
    }

    void ChElementShellTri_3::ComputeMmatrixGlobal(ChMatrix<>& M)
    {
        getElementMR(M, 1.0, 0.0);
    }


    void ChElementShellTri_3::SetNodes(std::shared_ptr<ChNodeFEAxyz> node0, std::shared_ptr<ChNodeFEAxyz> node1, std::shared_ptr<ChNodeFEAxyz> node2)
    {
        all_nodes[0] = node0;
        all_nodes[1] = node1;
        all_nodes[2] = node2;
    }

    void ChElementShellTri_3::SetNodes(std::shared_ptr<ChNodeFEAxyz> node0, std::shared_ptr<ChNodeFEAxyz> node1, std::shared_ptr<ChNodeFEAxyz> node2, std::shared_ptr<ChNodeFEAxyz> node3, std::shared_ptr<ChNodeFEAxyz> node4, std::shared_ptr<ChNodeFEAxyz> node5)
    {
        SetNodes(node0, node1, node2);
        SetNeighbouringNodes(node3, node4, node5);
    }

    void ChElementShellTri_3::SetNeighbouringNodes(std::shared_ptr<ChNodeFEAxyz> node3, std::shared_ptr<ChNodeFEAxyz> node4, std::shared_ptr<ChNodeFEAxyz> node5)
    {
        all_nodes[3] = node3;
        all_nodes[4] = node4;
        all_nodes[5] = node5;
    }

}  // END_OF_NAMESPACE____
}  // END_OF_NAMESPACE____
