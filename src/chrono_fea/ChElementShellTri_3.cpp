#include "chrono_fea/ChElementShellTri_3.h"
#include <bitset>

namespace chrono {
namespace fea {

    double ChElementShellTri_3::GetArea0() const
    {
        ChVector<double> vect_temp;
        vect_temp.Cross(GetEdgeVector(0, 0), GetEdgeVector(1, 0));
        return vect_temp.Length() / 2;
    }

    double ChElementShellTri_3::GetArea() const
    {
        ChVector<double> vect_temp;
        vect_temp.Cross(GetEdgeVector(0, 0), GetEdgeVector(1, 0));
        return vect_temp.Length() / 2;
    }

    void ChElementShellTri_3::GetEdgeVector0(ChVector<double>& vect, int edge_sel, int elem_sel) const
    {
        vect.Sub(all_nodes[GetNodeOfEdge(edge_sel, 1, elem_sel)]->GetX0(), all_nodes[GetNodeOfEdge(edge_sel, 0, elem_sel)]->GetX0());
    }

    ChVector<double> ChElementShellTri_3::GetEdgeVector0(int edge_sel, int elem_sel) const
    {
        ChVector<double> vect;
        vect.Sub(all_nodes[GetNodeOfEdge(edge_sel, 1, elem_sel)]->GetX0(), all_nodes[GetNodeOfEdge(edge_sel, 0, elem_sel)]->GetX0());
        return vect;
    }

    void ChElementShellTri_3::GetEdgeVector(ChVector<double>& vect, int edge_sel, int elem_sel) const
    {
        vect.Sub(all_nodes[GetNodeOfEdge(edge_sel, 1, elem_sel)]->GetPos(), all_nodes[GetNodeOfEdge(edge_sel, 0, elem_sel)]->GetPos());
    }

    ChVector<double> ChElementShellTri_3::GetEdgeVector(int edge_sel, int elem_sel) const
    {
        ChVector<double> vect;
        vect.Sub(all_nodes[GetNodeOfEdge(edge_sel, 1, elem_sel)]->GetPos(), all_nodes[GetNodeOfEdge(edge_sel, 0, elem_sel)]->GetPos());
        return vect;
    }

    double ChElementShellTri_3::GetEdgeVersor0(ChVector<double>& vers, int edge_sel, int elem_sel) const
    {
        GetEdgeVector0(vers, edge_sel, elem_sel);
        double norm = vers.Length();
        vers.Scale(1 / norm);
        return norm;
    }

    int ChElementShellTri_3::GetNodeOfEdge(int edge_sel, int node_sel, int elem_sel)
    {
        return edge_num[elem_sel][edge_sel + node_sel];
    }

    double ChElementShellTri_3::GetHeight(int edge_sel) const
    {
        return 2 * GetArea() / GetEdgeVector(edge_sel, 0).Length();
    }

    double ChElementShellTri_3::GetNeighbourHeight(int edge_sel) const
    {
        return 2 * neighbouring_elements[edge_sel]->GetArea() / GetEdgeVector(edge_sel, 0).Length();
    }


    void ChElementShellTri_3::initializeElement()
    {
        // store edge versors
        for (auto edge_sel = 0; edge_sel < 3; ++edge_sel)
        {
            edge_length0[edge_sel] = GetEdgeVersor0(edge_versors0[edge_sel], edge_sel, 0);
        }

        // compute area of main element
        element_area0 = GetArea0();


        // normal in-plane versors of the main element M (glob ref)
        // project the triangle on the XY plane
        // compute edge versors, rotate -90° around Z (outward normal)
        for (auto edge_sel = 0; edge_sel < 3; ++edge_sel)
        {
            double norm_temp = edge_versors0[edge_sel](0) * edge_versors0[edge_sel](0) + edge_versors0[edge_sel](1) * edge_versors0[edge_sel](1);
            ip_normal_mod[edge_sel](0) = edge_versors0[edge_sel](1) / norm_temp * edge_length0[edge_sel] / 2 / element_area0;
            ip_normal_mod[edge_sel](1) = -edge_versors0[edge_sel](0) / norm_temp * edge_length0[edge_sel] / 2 / element_area0;
        }

        // compute scale factors
        for (auto edge_sel = 0; edge_sel < 3; ++edge_sel)
        {
            for (auto col_sel = 0; col_sel < 3; ++col_sel)
            {
                J_source(edge_sel, col_sel) = edge_length0[col_sel] * ChVector<double>::Dot(edge_versors0[col_sel], edge_versors0[edge_sel]) / 2 / element_area0;
                if (neighbouring_elements[edge_sel])
                {
                    J_source(edge_sel, col_sel + 3) = -ChVector<double>::Dot(GetEdgeVector0(col_sel, edge_sel+1), edge_versors0[edge_sel]) / 2 / neighbouring_elements[edge_sel]->GetArea0();
                }
            }
        }

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
            edge_bc[2] = CLAMPED;

        if (all_nodes[1]->GetFixed() && all_nodes[2]->GetFixed())
            edge_bc[0] = CLAMPED;

        if (all_nodes[2]->GetFixed() && all_nodes[0]->GetFixed())
            edge_bc[1] = CLAMPED;
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

    void ChElementShellTri_3::Update()
    {
        /* ************** Compute the curvature matrix ************** */
        // compute r_i^M
        std::array<double, 3> stiff_ratio;
        for (auto edge_sel = 0; edge_sel < 3; ++edge_sel)
        {
            if (neighbouring_elements[edge_sel])
                stiff_ratio[edge_sel] = 1 / (1 + (m_material->Get_E() * pow(thickness, 3) * GetHeight(edge_sel)) / (neighbouring_elements[edge_sel]->m_material->Get_E() * pow(neighbouring_elements[edge_sel]->thickness, 3) * GetNeighbourHeight(edge_sel)));
            else if (edge_bc[edge_sel] == CLAMPED) // only clamped edge BC supported
                stiff_ratio[edge_sel] = 1;
            else
                stiff_ratio[edge_sel] = 0;
        }

        // compute normal to plane
        std::array<ChVector<double>, 4> t;
        for (auto elem_sel = 0; elem_sel < 4; ++elem_sel)
        {
            if (all_nodes[elem_sel + 2]) // check if neighbour exists before trying to compute area
            {
                t[elem_sel].Cross(GetEdgeVector(0, elem_sel), GetEdgeVector(1, elem_sel));
                t[elem_sel].Normalize();
            }
        }

        // compute J
        ChMatrixNM<double, 3, 18> B_bend;
        for (auto edge_sel = 0; edge_sel < 3; ++edge_sel)
        {
            for (auto col_sel = 0; col_sel < 3; ++col_sel)
            {
                B_bend.PasteSumVector(ip_normal_mod[edge_sel] * (stiff_ratio[edge_sel] * t[0](col_sel) * J_source(edge_sel, col_sel)), 0, col_sel * 3);
                B_bend.PasteSumVector(ip_normal_mod[edge_sel] * (stiff_ratio[edge_sel] * t[0](col_sel) * J_source(edge_sel, col_sel)), 0, col_sel * 3 + 1);
                B_bend.PasteSumVector(ip_normal_mod[edge_sel] * (stiff_ratio[edge_sel] * t[0](col_sel) * J_source(edge_sel, col_sel)), 0, col_sel * 3 + 2);

                if (all_nodes[edge_sel + 3]) // check if neighbour exists before trying to compute area
                {
                    B_bend.PasteSumVector(ip_normal_mod[edge_sel] * (stiff_ratio[edge_sel] * t[edge_sel + 1](col_sel) * J_source(edge_sel, col_sel)), 0, col_sel * 3);
                    B_bend.PasteSumVector(ip_normal_mod[edge_sel] * (stiff_ratio[edge_sel] * t[edge_sel + 1](col_sel) * J_source(edge_sel, col_sel)), 0, col_sel * 3 + 1);
                    B_bend.PasteSumVector(ip_normal_mod[edge_sel] * (stiff_ratio[edge_sel] * t[edge_sel + 1](col_sel) * J_source(edge_sel, col_sel)), 0, col_sel * 3 + 2);
                }
            }
        }

        /* ***************************************************** */


        // Update bending stiffness matrix
        ChMatrixNM<double, 18, 18> K_bend;
        ChMatrixNM<double, 18, 3> mat_temp;
        ChMatrixNM<double, 3, 3> temp = m_material->GetConsitutiveMatrixBending();
        mat_temp.MatrTMultiply(B_bend, temp);
        K_bend.MatrMultiply(mat_temp, B_bend);

        stiffness_matrix = K_bend;
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
        initializeElement();

        // Inform the system about which nodes take part to the computation
        //TODO: do I have to put also neighbouring nodes?
        std::vector<ChVariables*> vars;
        for (auto node_sel :all_nodes)
        {
            if (node_sel)
                vars.push_back(&node_sel->Variables());
        }

        Kmatr.SetVariables(vars);
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
