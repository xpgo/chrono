#include "chrono_fea/ChElementShellTri_3.h"
#include <bitset>

namespace chrono {
namespace fea {

    double ChElementShellTri_3::GetArea0()
    {
        if (element_area0<0)
        {
            ChVector<double> vect_temp;
            vect_temp.Cross(GetEdgeVector0(0, -1), GetEdgeVector0(1, -1));
            element_area0 = vect_temp.Length() / 2;
        }

        std::cout << element_area0 << std::endl;
        
        return element_area0;
    }

    double ChElementShellTri_3::GetArea() const
    {
        ChVector<double> vect_temp;
        vect_temp.Cross(GetEdgeVector(0, -1), GetEdgeVector(1, -1));
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
        return edge_num[elem_sel+1][edge_sel + node_sel];
    }

    double ChElementShellTri_3::GetHeight(int edge_sel) const
    {
        return 2 * GetArea() / GetEdgeVector(edge_sel, -1).Length();
    }

    double ChElementShellTri_3::GetNeighbourHeight(int edge_sel) const
    {
        return 2 * neighbouring_elements[edge_sel]->GetArea() / GetEdgeVector(edge_sel, -1).Length();
    }


    void ChElementShellTri_3::initializeElement()
    {

        ChVector<double> vectX;
        ChVector<double> vectY;
        ChVector<double> vectZ;




        // store edge versors
        for (auto edge_sel = 0; edge_sel < 3; ++edge_sel)
        {
            edge_length0[edge_sel] = GetEdgeVersor0(edge_versors0[edge_sel], edge_sel, -1);
        }

        
        // Shape function derivatives
        ChMatrixNM<double, 2, 3> shape_fun_der;
        for (auto edge_sel = 0; edge_sel < 3; ++edge_sel)
        {
            // TODO: this method is WRONG
            // normal in-plane versors of the main element M (glob ref)
            // project the triangle on the XY plane
            // compute edge versors, rotate -90° around Z (outward normal)
            double norm_temp = edge_versors0[edge_sel](0) * edge_versors0[edge_sel](0) + edge_versors0[edge_sel](1) * edge_versors0[edge_sel](1);
            shape_fun_der(0, edge_sel) = edge_versors0[edge_sel](1) / norm_temp * edge_length0[edge_sel] / 2 / element_area0;
            shape_fun_der(1, edge_sel) = -edge_versors0[edge_sel](0) / norm_temp * edge_length0[edge_sel] / 2 / element_area0;

            L_block(0, edge_sel) = pow(shape_fun_der(0, edge_sel), 2);
            L_block(1, edge_sel) = pow(shape_fun_der(1, edge_sel), 2);
            L_block(2, edge_sel) = -2* shape_fun_der(0, edge_sel)*shape_fun_der(1, edge_sel);
        }

        // derivative of phi (considered constant)
        ChMatrixNM<double, 3, 2> phi_der;
        for (auto der_sel = 0; der_sel < 2; ++der_sel)
        {
            for (auto node_sel = 0; node_sel < 3; ++node_sel)
            {
                phi_der.PasteSumVector(all_nodes[node_sel]->GetPos() * shape_fun_der(der_sel, node_sel), 0, der_sel);
            }
        }

        ChMatrixNM<double, 3, 18> Bmembr;
        for (auto row_sel = 0; row_sel<3; ++row_sel)
        {
            ChMatrixNM<double, 3, 1> mat_temp{ phi_der };
            for (auto col_sel = 0; col_sel<3; ++col_sel)
            {
            }
        }
        
        

        // compute main projectors
        for (auto edge_sel = 0; edge_sel < 3; ++edge_sel)
        {
            for (auto node_sel = 0; node_sel < 3; ++node_sel)
            {
                main_projectors(edge_sel, node_sel) = edge_length0[node_sel] * ChVector<double>::Dot(edge_versors0[node_sel], edge_versors0[edge_sel]) / 2 / element_area0;
            }
        }

        // compute main projectors
        for (auto edge_sel = 0; edge_sel < 3; ++edge_sel)
        {
            for (auto node_sel = 0; node_sel < 3; ++node_sel)
            {
                if (neighbouring_elements[edge_sel])
                {
                    neigh_projectors(edge_sel, node_sel) = -ChVector<double>::Dot(GetEdgeVector0(node_sel, edge_sel), edge_versors0[edge_sel]) / 2 / neighbouring_elements[edge_sel]->GetArea0();
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

    void ChElementShellTri_3::GetElementData(const ChVector<>& P1_in, const ChVector<>& P2_in, const ChVector<>& P3_in,
                                             ChMatrix<>* gradient_local, ChVector<>* z_versor, double area, ChMatrix<>* gradient_shape_function, ChMatrix<>* edge_normal_vers, std::array<double,3>* edge_length, ChMatrix33<>* rotGL, ChMatrix<>* c_proj, ChMatrix<>* gradient_side_n)
    {
        //ChVector<double> z_versor;
        //std::array<double, 3> edge_length;
        //ChMatrix33<double> rotGL;
        //ChMatrixNM<double, 2, 3> edge_normal_vers;
        //ChMatrixNM<double, 3, 2> gradient_shape_function;
        //ChMatrix33<double> c_proj;
        //ChMatrixNM<double,3,2> gradient_local;



        // edge vectors
        ChVector<double> e_vect1;
        ChVector<double> e_vect2;
        ChVector<double> e_vect3;

        e_vect1.Sub(P1_in, P2_in);
        e_vect3.Sub(P1_in, P3_in);
        e_vect3.Sub(P2_in, P1_in);


        (*edge_length)[0] = e_vect1.Length();
        (*edge_length)[1] = e_vect2.Length();
        (*edge_length)[2] = e_vect3.Length();

        // local Cartesian axis
        ChVector<double> x_versor;
        ChVector<double> y_versor;
        
        z_versor->Cross(e_vect3, -e_vect2);
        auto lambda = z_versor->Length();
        area = lambda / 2; //TODO: strange fact...
        z_versor->Dot(1 / lambda);
        x_versor = e_vect3;
        x_versor.Normalize();
        y_versor.Cross(*z_versor, x_versor);

        // rotation matrix from Local coordinates to Global coordinates
        rotGL->PasteVector(x_versor, 0, 0);
        rotGL->PasteVector(y_versor, 0, 1);
        rotGL->PasteVector(*z_versor, 0, 2);

        ChMatrix33<double> rotLG;
        rotLG.FastInvert(rotGL); //TODO: can be done in-place on rotGL?

        ChMatrix33<double> nodes_local;
        //nodes_local.PasteVector(rotLG.Matr_x_Vect(P1_in - P1_in) , 0, 0); 
        nodes_local.PasteVector(rotLG.Matr_x_Vect(P2_in - P1_in), 0, 0);
        nodes_local.PasteVector(rotLG.Matr_x_Vect(P3_in - P1_in), 0, 0);

        // normal-to-edge vectors
        ChVector<double> n_vect1;
        ChVector<double> n_vect2;
        ChVector<double> n_vect3;

        ChMatrixNM<double, 2, 3> edge_normal_vect;
        edge_normal_vect(0, 0) = + nodes_local(1, 2) - nodes_local(1, 1);
        edge_normal_vect(0, 1) = + nodes_local(1, 0) - nodes_local(1, 2);
        edge_normal_vect(0, 2) = + nodes_local(1, 1) - nodes_local(1, 0);
        edge_normal_vect(1, 0) = - nodes_local(0, 2) + nodes_local(0, 1);
        edge_normal_vect(1, 1) = - nodes_local(0, 0) + nodes_local(0, 2);
        edge_normal_vect(1, 2) = - nodes_local(0, 1) + nodes_local(0, 0);

        // gradient shape function
        gradient_shape_function->CopyFromMatrixT(edge_normal_vect);
        *gradient_shape_function *= -1/2/area;


        // gradients (local frame)
        ChMatrix33<double> nodes_glob;
        nodes_glob.PasteVector(P1_in, 0, 0);
        nodes_glob.PasteVector(P2_in, 0, 1);
        nodes_glob.PasteVector(P3_in, 0, 2);
        gradient_local->MatrMultiply(nodes_glob, *gradient_shape_function);

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

        // gradients (side frame)
        ChMatrixNM<double, 2, 1> edge_normal_vers_temp;
        edge_normal_vers_temp(0, 0) = -(*edge_normal_vers)(0, 2);
        edge_normal_vers_temp(1, 0) = -(*edge_normal_vers)(1, 2);
        gradient_side_n->MatrMultiply(*gradient_local, edge_normal_vers_temp);

        // C projections
        ChMatrix33<double> c_temp1;
        c_temp1.PasteVector(e_vect1, 0, 0);
        c_temp1.PasteVector(e_vect2, 0, 1);
        c_temp1.PasteVector(e_vect3, 0, 2);

        ChMatrix33<double> c_temp2;
        c_temp2.PasteVector(e_vect1/ (*edge_length)[0], 0, 0);
        c_temp2.PasteVector(e_vect2/ (*edge_length)[1], 0, 1);
        c_temp2.PasteVector(e_vect3/ (*edge_length)[2], 0, 2);

        c_proj->MatrTMultiply(c_temp1, c_temp2);
        *c_proj *= 1 / 2 / area;
        

    }

    void ChElementShellTri_3::Update()
    {
        ChMatrixNM<double,18,1> lambda_gamma;
        ChMatrixNM<double,3,1> L_main_undeformed;
        ChMatrixNM<double, 3, 18> Bb;
        ChMatrixNM<double, 3, 18> Bm;

        GetElementData(all_nodes[0]->GetPos(), all_nodes[1]->GetPos(), all_nodes[2]->GetPos(),
                       &gradient_local, &z_versor, element_area, &gradient_shape_function, &edge_normal_vers, &edge_length, &rotGL, &c_proj, &gradient_side_n);
        update_counter++;
        

        for (auto edge_sel = 0; edge_sel<3; edge_sel++)
        {
            // compute lamba*gamma
            lambda_gamma.PasteVector(z_versor*c_proj0(0,edge_sel), 0, 0);
            lambda_gamma.PasteVector(z_versor*c_proj0(1,edge_sel), 3, 0);
            lambda_gamma.PasteVector(z_versor*c_proj0(2,edge_sel), 6, 0);

            for (auto node_sel = 0; node_sel<3; node_sel++)
            {
                if (neighbouring_elements[node_sel])
                {
                    lambda_gamma.PasteSumVector(neighbouring_elements[node_sel]->z_versor*neighbouring_elements[node_sel]->c_proj0(node_sel, 2),node_sel+3,0);
                }
            }

            lambda_gamma.MatrTranspose();

            // compute L block
            L_main_undeformed(0, 0) = pow(gradient_shape_function0(edge_sel, 0), 2);
            L_main_undeformed(1, 0) = pow(gradient_shape_function0(edge_sel, 1), 2);
            L_main_undeformed(2, 0) = 2*gradient_shape_function0(edge_sel, 0)*gradient_shape_function0(edge_sel, 1);

            // relative stiffness
            // TODO: relate rM to Young mod and thickness
            double riM = 0.5;

            // curvature matrix
            
            ChMatrixNM<double, 3, 18> Bb_temp;
            Bb_temp.MatrMultiply(L_main_undeformed*(riM / edge_length0[edge_sel]), lambda_gamma);
            Bb.PasteSumMatrix(&Bb_temp,0,0);

            ChMatrixNM<double, 3, 3> Bm_temp;
            Bb_temp.PasteVector(gradient_local.ClipVector(0, 0)*gradient_shape_function0(edge_sel, 0), 0, 0);
            Bb_temp.PasteVector(gradient_local.ClipVector(0, 1)*gradient_shape_function0(edge_sel, 1), 0, 1);
            Bb_temp.PasteVector(gradient_local.ClipVector(0, 1)*gradient_shape_function0(edge_sel, 0) +
                                gradient_local.ClipVector(0, 0)*gradient_shape_function0(edge_sel, 1), 0, 2);
            Bb.PasteTranspMatrix(&Bb_temp, 0, edge_sel * 3);
        }

        
        Bb *= 4 * element_area0;
        
        // constitutive matrix //TODO: to be moved elsewhere
        double nu = 0.5;
        double YoungMod = 210e9;
        ChMatrixNM<double, 3, 3> constitutive_matrix;
        constitutive_matrix(0, 0) = 1;
        constitutive_matrix(0, 1) = nu;
        constitutive_matrix(1, 0) = nu;
        constitutive_matrix(1, 1) = 1;
        constitutive_matrix(2, 2) = (1 - nu)/2;
        constitutive_matrix *= YoungMod / (1 - nu*nu);


        // stiffness matrix
        ChMatrixNM<double, 3, 18> temp;
        temp.MatrMultiply(constitutive_matrix, Bm);
        temp *= thickness0;
        stiffness_matrix.MatrTMultiply(Bm, temp);
        temp.MatrMultiply(constitutive_matrix, Bb);
        temp *= pow(thickness0,3)/12;
        ChMatrixNM<double, 18, 18> stiffness_matrix_membr;
        stiffness_matrix_membr.MatrTMultiply(Bm, temp);
        stiffness_matrix += stiffness_matrix_membr;



    }


    void ChElementShellTri_3::computeLt(ChMatrix<double>& mat_temp, int elem_sel, std::array<ChVector<double>, 4>& t, double scale)
    {
        for (auto row_sel = 0; row_sel<3; ++row_sel)
        {
            for (auto col_sel = 0; col_sel<3; ++col_sel)
            {
                mat_temp(row_sel, col_sel) = scale*L_block(row_sel, 0)*t[elem_sel+1](col_sel);
            }
        }
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
        GetElementData(all_nodes[0]->GetX0(), all_nodes[1]->GetX0(), all_nodes[2]->GetX0(),
                       &gradient_local0, &z_versor0, element_area0, &gradient_shape_function0, &edge_normal_vers0, &edge_length0, &rotGL0, &c_proj0, &gradient_side_n0);

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
