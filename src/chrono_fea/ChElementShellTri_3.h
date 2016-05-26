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
// File authors: Dario Mangoni, Alessandro Tasora

#ifndef CHELEMENTSHELLTRI_H
#define CHELEMENTSHELLTRI_H

#include "chrono_fea/ChNodeFEAxyz.h"
#include "chrono_fea/ChElementShell.h"
#include "ChMesh.h"

namespace chrono {
namespace fea {

/// @addtogroup fea_elements
/// @{

/// Material definition.
    class ChApiFea ChMaterialShellTri_3 {
    public:
        /// Construct an isotropic material.
        ChMaterialShellTri_3(double E,    ///< Young's modulus
            double nu,   ///< Poisson ratio
            double rho
        ) :
            m_rho(rho),
            m_E(E),
            m_nu(nu)
        {
            UpdateConsitutiveMatrices();
        }

        /// Return the elasticity moduli
        double Get_E() const { return m_E; }
        void Set_E(double E_in) { m_E = E_in; UpdateConsitutiveMatrices(); }

        /// Return the Poisson ratio
        double Get_nu() const { return m_nu; }
        void Set_nu(double nu_in) { m_nu = nu_in; UpdateConsitutiveMatrices();
        }

        /// Return the density
        double Get_rho() const { return m_rho; }
        void Set_rho(double rho_in) { m_rho = rho_in; UpdateConsitutiveMatrices();
        }

        void UpdateConsitutiveMatrices()
        {
            updateDbend();
            updateDmembr();
        }

        ChMatrixNM<double, 3, 3>& GetConsitutiveMatrixBending() { return Dbend; };
        ChMatrixNM<double, 3, 3>& GetConsitutiveMatrixMembrane() { return Dmembr; };


        /// The FE code will evaluate this function to compute 
        /// u,v stresses/torques given the u,v strains/curvatures.
        /// You can inherit a more sophisticated material that override this (ex. for
        /// orthotropic materials, etc.)

        //virtual void ComputeStress(ChVector<>& n_u,
        //    ChVector<>& n_v,
        //    ChVector<>& m_u,
        //    ChVector<>& m_v,
        //    const ChVector<>& eps_u,
        //    const ChVector<>& eps_v,
        //    const ChVector<>& kur_u,
        //    const ChVector<>& kur_v);

    protected:
        void updateDbend()
        {
            double Dbend_multiplier = m_E / 12 / (1 - m_nu*m_nu);
            Dbend(0, 0) = Dbend_multiplier;
            Dbend(0, 1) = m_nu*Dbend_multiplier;
            Dbend(1, 0) = m_nu*Dbend_multiplier;
            Dbend(1, 1) = Dbend_multiplier;
            Dbend(2, 2) = Dbend_multiplier*(1-m_nu)/2;
        }

        void updateDmembr()
        {
            double Dmembr_multiplier = m_E / (1 - m_nu*m_nu);
            Dmembr(0, 0) = Dmembr_multiplier;
            Dmembr(0, 1) = m_nu*Dmembr_multiplier;
            Dmembr(1, 0) = m_nu*Dmembr_multiplier;
            Dmembr(1, 1) = Dmembr_multiplier;
            Dmembr(2, 2) = Dmembr_multiplier*(1 - m_nu) / 2;
        }

    private:

        ChMatrixNM<double, 3, 3> Dbend;
        ChMatrixNM<double, 3, 3> Dmembr;
        double m_rho = 7850;    ///< density
        double m_E = 210e9; ///< elasticity moduli
        double m_nu = 0.3;  ///< Poisson ratio
    };





/// Class for FEA elements of Onate shell type.
/// This element has a linear displacement field.
class ChApiFea ChElementShellTri_3 : public ChElementShell {
private:
    std::vector<std::shared_ptr<ChNodeFEAxyz>> main_nodes; ///< element nodes


    double edge_length[3];
    double element_area = 0;
    double thickness = 0;
    ChMatrixNM<double, 18, 18> stiffness_matrix;

    size_t updates_count = 0;

    /// edge versors;
    /// s_loc[1] is of edge 1, that goes from 2 to 3
    /// s_loc[2] is of edge 2, that goes from 3 to 1
    /// s_loc[3] is of edge 3, that goes from 1 to 2
    std::vector<ChMatrixNM<double, 2, 1>> s_loc;

    ChMatrix33<double> rotGL; ///< rotation matrix from Local to Global frame

    enum boundary_conditions
    {
        DEFAULT = 0,
        CLAMPED = 1,
        SUPPORTED = 2,
        SYMMETRIC = 3
    } edge_bc[3] = { DEFAULT, DEFAULT, DEFAULT };

    std::vector<std::shared_ptr<ChElementShellTri_3>> neighbouring_elements; ///< neighbour elements

    /// column j-th refers to neighbour element attached to edge 'j'
    /// row i-th refers to the i-th node of the neighbouring element
    /// the array tells that the i-th node of the j-th neighbour
    /// is the neighbour_nodes[i][j]-th node of the current element
    /// if neighbour_nodes[i][j] < 3 the node belongs to the current element
    /// if neighbour_nodes[i][j] >= 3 the node actually belongs to the j-th neighbour
    /// so that the neighbour attached to edge1 has node 4
    /// the neighbour attached to edge2 has node 5
    /// the neighbour attached to edge3 has node 6
    int neighbour_nodes[3][3] = {{-1, -1, -1},{-1, -1 ,-1},{ -1, -1 ,-1 } };

    /// column j-th refers to neighbour element attached to edge 'j'
    /// if neighbour_node_not_shared[i][j] < 3 the node belongs to the current element
    /// if neighbour_node_not_shared[i][j] >= 3 the node actually belongs to the j-th neighbour
    /// so that the neighbour attached to edge1 has node 4
    /// the neighbour attached to edge2 has node 5
    /// the neighbour attached to edge3 has node 6
    int neighbour_node_not_shared[3] = { -1,-1,-1 };
    std::shared_ptr<ChMaterialShellTri_3> m_material;

    ChMatrix33<double> shape_function;
    ChMatrixNM<double, 9, 1> minus_pos_init; ///< stores the initial position of the three main nodes, sign changed

    int countNeighbours() const
    {
        int neighbours_count = 0;;
        for (auto neigh_node_sel = 0; neigh_node_sel < 3; neigh_node_sel++)
        {
            if (neighbour_node_not_shared[neigh_node_sel] != -1)
                neighbours_count++;
        }
        return neighbours_count;
    }

    void updateGeometry() {

        updates_count++;

        // evaluate x' (x local axis) in global coordinates
        ChVector<double> x_loc;
        x_loc = main_nodes[1]->GetPos() - main_nodes[0]->GetPos();
        x_loc.Normalize();
        rotGL.PasteVector(x_loc, 0, 0);
        // evaluate y' (y local axis) in global coordinates
        ChVector<double> y_loc;
        y_loc = main_nodes[2]->GetPos() - main_nodes[0]->GetPos();
        double proj_y_on_x = y_loc.Dot(x_loc);
        y_loc.Sub(y_loc, proj_y_on_x*x_loc);
        y_loc.Normalize();
        rotGL.PasteVector(y_loc, 0, 1);
        // evaluate z' (z local axis) in global coordinates
        ChVector<double> z_loc;
        z_loc.Cross(x_loc, y_loc);
        rotGL.PasteVector(z_loc, 0, 2);

        ChMatrix33<double> rotLG;
        rotGL.FastInvert(&rotLG);
        // local position of nodes
        // - takes the global position vector of the three main nodes
        // - rotates into the local frame
        // - writes the x and y components into temp_shape
        // - invert temp_shape to get the shape function
        ChMatrix33<double> temp_shape;
        std::vector<ChMatrixNM<double, 2, 1>> nodes_pos_loc;
        nodes_pos_loc.resize(3);
        ChVector<double> temp_vect;
        for (auto node_sel = 0; node_sel < 3; node_sel++)
        {
            temp_vect = rotLG*main_nodes[node_sel]->GetPos();
            temp_shape(node_sel, 0) = 1;
            temp_shape(node_sel, 1) = temp_vect(0);
            temp_shape(node_sel, 2) = temp_vect(1);
            nodes_pos_loc[node_sel](0, 0) = temp_vect(0);
            nodes_pos_loc[node_sel](1, 0) = temp_vect(1);
        }

        element_area = temp_shape.FastInvert(&shape_function) / 2;

        // compute the directions of the edges (not normalized yet)
        s_loc[0] = nodes_pos_loc[2] - nodes_pos_loc[1]; // edge1 goes from 2 to 3
        s_loc[1] = nodes_pos_loc[0] - nodes_pos_loc[2]; // edge2 goes from 3 to 1
        s_loc[2] = nodes_pos_loc[1] - nodes_pos_loc[0]; // edge3 goes from 1 to 2

        // normalize the 's' vectors and build 't' versors
        for (auto edge_sel = 0; edge_sel < 3; edge_sel++)
        {
            edge_length[edge_sel] = s_loc[edge_sel].NormTwo();
            s_loc[edge_sel].MatrDivScale(edge_length[edge_sel]);
        }


    }


    void updateStructural()
    {

        for (auto neigh_sel = 0; neigh_sel < 3; neigh_sel++)
        {
            if (neighbouring_elements[neigh_sel].get() != nullptr && neighbouring_elements[neigh_sel]->updates_count < updates_count)
                neighbouring_elements[neigh_sel]->updateGeometry();
        }

        ///////////////////// Bending stiffness //////////////////////////
        ChMatrixNM<double, 3, 12> Bll_bend; ///< Curvature Matrix (that takes local displacements and returns local curvatures)
        ChMatrixNM<double, 2, 12> J;

        for (auto edge_sel = 0; edge_sel < 3; edge_sel++)
        {
            J.FillElem(0);

            // theta_s part
            if (!(edge_bc[edge_sel] == CLAMPED || edge_bc[edge_sel] == SYMMETRIC))
            {
                if (neighbouring_elements[edge_sel].get()!=nullptr)
                {
                    // fixed part of theta_s
                    for (auto col_sel = 0; col_sel<3; col_sel++)
                        J(0, col_sel) = 0.5 * (-s_loc[edge_sel](1) * shape_function(1, col_sel) + s_loc[edge_sel](0) * shape_function(2, col_sel));

                    // moving part of theta_s
                    for (auto col_sel = 0; col_sel<3; col_sel++)
                        J(0, (edge_sel+1)*3+col_sel) = -0.5 * (-neighbouring_elements[edge_sel]->s_loc[neighbour_node_not_shared[edge_sel]](1) * neighbouring_elements[edge_sel]->shape_function(1, col_sel)
                                                               +neighbouring_elements[edge_sel]->s_loc[neighbour_node_not_shared[edge_sel]](0) * neighbouring_elements[edge_sel]->shape_function(2, col_sel) );
                }
                else
                {
                    // fixed part of theta_s
                    for (auto col_sel = 0; col_sel<3; col_sel++)
                        J(0, col_sel) = 1 * (-s_loc[edge_sel](1) * shape_function(1, col_sel) + s_loc[edge_sel](0) * shape_function(2, col_sel));
                }
            }
                

            // theta_t part
            if (!( edge_bc[edge_sel]==CLAMPED || edge_bc[edge_sel] == SUPPORTED ))
            {
                switch (edge_sel)
                {
                    case 0:
                        J(1, 1) = +1.0 / edge_length[edge_sel];
                        J(1, 2) = -1.0 / edge_length[edge_sel];
                        break;
                    case 1:
                        J(1, 0) = -1.0 / edge_length[edge_sel];
                        J(1, 2) = +1.0 / edge_length[edge_sel];
                        break;
                    case 2:
                        J(1, 0) = +1.0 / edge_length[edge_sel];
                        J(1, 1) = -1.0 / edge_length[edge_sel];
                        break;
                    default:
                        assert(0);
                }
                
            }

            // This temporary matrix includes T * [0, -1; 1, 0] * rotLS2D
            ChMatrixNM<double, 3, 2> rot_temp;
            rot_temp(0, 0) = -s_loc[edge_sel](1)*s_loc[edge_sel](1);
            rot_temp(0, 1) = -s_loc[edge_sel](0)*s_loc[edge_sel](1);
            rot_temp(1, 0) = -s_loc[edge_sel](0)*s_loc[edge_sel](0);
            rot_temp(1, 1) = -rot_temp(0, 1);
            rot_temp(2, 0) = 2 * rot_temp(1, 1);
            rot_temp(2, 1) = rot_temp(0, 0) - rot_temp(1, 0);

            // Update the Curvature Matrix 'Bll'
            ChMatrixNM<double, 3, 12> Bll_temp;
            Bll_temp.MatrMultiply(rot_temp, J);
            Bll_temp.MatrScale(edge_length[edge_sel]);
            Bll_bend.PasteSumMatrix(&Bll_temp,0,0);

        }

        Bll_bend.MatrDivScale(element_area);

        ChMatrixNM<double, 18, 12> rotLGw_transp;
        // first: fix the main element rotation part
        rotLGw_transp.PasteClippedMatrix(&rotGL, 0, 2, 3, 1, 0, 0);
        rotLGw_transp.PasteClippedMatrix(&rotGL, 0, 2, 3, 1, 3, 1);
        rotLGw_transp.PasteClippedMatrix(&rotGL, 0, 2, 3, 1, 6, 2);

        for (auto neigh_elem_sel = 0; neigh_elem_sel < 3; neigh_elem_sel++)
        {
            if (neighbouring_elements[neigh_elem_sel].get() == nullptr)
                continue;

            for (auto neigh_node_sel = 0; neigh_node_sel < 3; neigh_node_sel++)
            {
                rotLGw_transp.PasteClippedMatrix(&(neighbouring_elements[neigh_elem_sel]->rotGL), 0, 2, 3, 1, 3 * neighbour_nodes[neigh_node_sel][neigh_elem_sel], 3 * (neigh_elem_sel + 1) + neigh_node_sel);
            }

        }

        ChMatrixNM<double, 3, 18> Blg_bend;
        // Update the Curvature Matrix that takes Global displacements and returns Local curvatures
        Blg_bend.MatrMultiplyT(Bll_bend, rotLGw_transp);

        // Update bending stiffness matrix
        ChMatrixNM<double, 18, 18> K_bend;
        ChMatrixNM<double, 18, 3> mat_temp;
        ChMatrixNM<double, 3, 3> test = m_material->GetConsitutiveMatrixBending();
        mat_temp.MatrTMultiply(Blg_bend, test);
        mat_temp.MatrScale(pow(thickness, 3)*element_area);
        K_bend.MatrMultiply(mat_temp, Blg_bend);

        ///////////////////// Membrane stiffness //////////////////////////
        ChMatrixNM<double, 18, 3> Blg_membr_transp;
        
        ChVector<double> temp = rotGL.ClipVector(0, 0);
        for (auto col_sel = 0; col_sel < 3; col_sel++)
        {
            Blg_membr_transp.PasteVector(temp*shape_function(1, col_sel), col_sel * 3, 0);
            Blg_membr_transp.PasteSumVector(temp*shape_function(2, col_sel), col_sel * 3, 2);
        }

        temp = rotGL.ClipVector(0, 1);
        for (auto col_sel = 0; col_sel < 3; col_sel++)
        {
            Blg_membr_transp.PasteVector(temp*shape_function(2, col_sel), col_sel * 3, 1);
            Blg_membr_transp.PasteSumVector(temp*shape_function(1, col_sel), col_sel * 3, 2);
        }

        ChMatrixNM<double, 18, 18> K_membr;
        mat_temp.MatrMultiply(Blg_membr_transp, m_material->GetConsitutiveMatrixMembrane());
        mat_temp.MatrScale(thickness*element_area);
        K_membr.MatrMultiplyT(mat_temp, Blg_membr_transp);


        // Compose stiffness matrix
        stiffness_matrix.MatrAdd(K_membr, K_bend);


    }

    void updateBC() //TODO: find other ways to implement boundary conditions
    {
        if (main_nodes[0]->GetFixed() && main_nodes[1]->GetFixed())
            edge_bc[2] = CLAMPED;

        if (main_nodes[1]->GetFixed() && main_nodes[2]->GetFixed())
            edge_bc[0] = CLAMPED;

        if (main_nodes[2]->GetFixed() && main_nodes[0]->GetFixed())
            edge_bc[1] = CLAMPED;
    }

    void updateK()
    {
        updateGeometry();
        updateStructural();
    }

    void updateElementMass(int node_sel)
    {
        mass = element_area * thickness * m_material->Get_rho();
        for (auto diag_sel = 0; diag_sel < 3; diag_sel++)
        {
            main_nodes[diag_sel]->SetMass(mass/3);
        }
    }

    void getNodeM(ChMatrix<double>& node_mass, int node_sel, double factor = 1.0)
    {
        if (node_mass.GetColumns() != 3 || node_mass.GetRows() != 3)
            node_mass.Resize(3, 3);

        updateElementMass(node_sel);
        node_mass.FillDiag(main_nodes[node_sel]->GetMass()*factor);

    }

    void getNodeK(ChMatrix<double>& node_damp, int node_sel, double factor = 1.0)
    {
        if (node_damp.GetColumns() != 3 || node_damp.GetRows() != 3)
            node_damp.Resize(3, 3);

        node_damp.FillDiag(1e-3);

    }

    void getElementMR(ChMatrix<>& H, double Mfactor, double Kfactor)
    {
         H.Resize(GetNdofs(), GetNdofs());

        // Fill the H matrix with damping and mass; K and M supposed diagonal-block
        ChMatrix33<double> node_mat;
        for (auto main_node_sel = 0; main_node_sel<3; main_node_sel++)
        {
            getNodeK(node_mat, main_node_sel, Kfactor);
            H.PasteSumMatrix(&node_mat, main_node_sel * 3, main_node_sel * 3);
            getNodeM(node_mat, main_node_sel, Mfactor);
            H.PasteSumMatrix(&node_mat, main_node_sel * 3, main_node_sel * 3);
        }

        int offset_diag = 0;
        for (auto diag_sel = 0; diag_sel<3; ++diag_sel)
        {
            if (neighbouring_elements[diag_sel].get() == nullptr)
            {
                offset_diag++;
                continue;
            }

            neighbouring_elements[diag_sel]->getNodeK(node_mat, neighbour_node_not_shared[diag_sel], Kfactor);
            H.PasteSumMatrix(&node_mat, (diag_sel + 3 - offset_diag) * 3, (diag_sel + 3 - offset_diag) * 3);
            neighbouring_elements[diag_sel]->getNodeM(node_mat, neighbour_node_not_shared[diag_sel], Mfactor);
            H.PasteSumMatrix(&node_mat, (diag_sel + 3 - offset_diag) * 3, (diag_sel + 3 - offset_diag) * 3);
        }

    }


public:
    ChElementShellTri_3() {
        main_nodes.resize(3);
        neighbouring_elements.resize(3);
        s_loc.resize(3);
    }

    /// Return the thickness
    double Get_thickness() const { return thickness; }
    void Set_thickness(double thickness_in) { thickness = thickness_in; }

    void SetNodes(std::shared_ptr<ChNodeFEAxyz> nodeA,
                  std::shared_ptr<ChNodeFEAxyz> nodeB,
                  std::shared_ptr<ChNodeFEAxyz> nodeC)
    {
        main_nodes[0] = nodeA;
        main_nodes[1] = nodeB;
        main_nodes[2] = nodeC;
    }

    void Update() override {
        updateGeometry();
        updateStructural();
    }

    void UpdateConnectivity(std::shared_ptr<ChMesh> mesh)
    {
        bool skip_this_element;
        int notshared_node;
        int neighbour_nodes_temp[3];
        int neighbours_found = 0;

        for (auto elem_sel = 0; elem_sel < mesh->GetNelements(); elem_sel++) // pick one element from the list
        {
            if (neighbours_found > 2)
                break;

            auto candidate = std::dynamic_pointer_cast<ChElementShellTri_3>(mesh->GetElement(elem_sel));
            skip_this_element = false;
            notshared_node = -1;
            neighbour_nodes_temp[0] = -1;
            neighbour_nodes_temp[1] = -1;
            neighbour_nodes_temp[2] = -1;

            for (auto neigh_node_sel = 0; neigh_node_sel < 3; neigh_node_sel++) // pick one node of the candidate neighbour
            {
                int main_node_sel;
                for (main_node_sel = 0; main_node_sel < 3; main_node_sel++) // pick one node of the current element
                {
                    if (candidate->GetNodeN(neigh_node_sel) == main_nodes[main_node_sel]) // if they match...
                    {
                        // store the information that the node 'neigh_node_sel' of the neighbour
                        // is actually the node 'main_node_sel' of the current element
                        neighbour_nodes_temp[neigh_node_sel] = main_node_sel;
                        break;
                    }
                } // main_node_sel loop

                if (main_node_sel==3) // the neighbour node hasn't been found within main_nodes
                {
                    if (notshared_node == -1) // ...give the candidate one more possibility: it might be exactly the non-shared node
                    {
                        notshared_node = neigh_node_sel;
                    }
                    else // ... but if the same candidate already throw away its possibility then skip it!
                    {
                        skip_this_element = true;
                        break; // --> it jumps into the neigh_node_sel loop
                    }
                }

            } // neigh_node_sel loop

            if (!skip_this_element && candidate.get()!=this)
            {
                // the candidate is actually a neighbour; but on which edge?
                // search the index of the main_node that is not shared with this element;
                // that would be the edge number on which the just-found neighbour is attached
                for (auto main_node_sel = 0; main_node_sel<3; ++main_node_sel)
                {
                    auto neigh_node_sel = 0;
                    for (; neigh_node_sel<3; ++neigh_node_sel)
                    {
                        if (neighbour_nodes_temp[neigh_node_sel] == main_node_sel)
                            break;
                    }

                    // if in 'neighbour_nodes' there is the point pointed to main_node_sel
                    // then the neigh_node_sel will be <3
                    if (neigh_node_sel == 3)
                    {
                        // the main_node_sel is not found in the 
                        neighbouring_elements[main_node_sel] = candidate;
                        neighbour_nodes_temp[notshared_node] = main_node_sel + 3;

                        neighbour_nodes[0][main_node_sel] = neighbour_nodes_temp[0];
                        neighbour_nodes[1][main_node_sel] = neighbour_nodes_temp[1];
                        neighbour_nodes[2][main_node_sel] = neighbour_nodes_temp[2];
                        neighbour_node_not_shared[main_node_sel] = notshared_node;

                        neighbours_found++;
                        break;
                    }
                        
                }

            }
        } // elem_sel loop

    }
    
    size_t GetUpdatesCount() const { return updates_count; }


    void SetBC(int BC_on_edge1, int BC_on_edge2, int BC_on_edge3)
    {
        edge_bc[0] = static_cast<boundary_conditions>(BC_on_edge1);
        edge_bc[1] = static_cast<boundary_conditions>(BC_on_edge2);
        edge_bc[2] = static_cast<boundary_conditions>(BC_on_edge3);
    }

    ChMatrixNM<double, 3, 3>& GetShapeFunction() { return shape_function; }


    

    void ComputeKRMmatricesGlobal(ChMatrix<>& H, double Kfactor, double Rfactor, double Mfactor) override {
        
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

    void ComputeInternalForces(ChMatrixDynamic<>& Fi) override {

        ChMatrixNM<double, 18, 1> disp_glob;
        ChMatrixNM<double, 18, 1> Fi_temp;


        updateK();


        // paste main_nodes position vectors
        for (auto node_sel = 0; node_sel<3; node_sel++)
        {
            disp_glob.PasteVector(main_nodes[node_sel]->GetPos(), node_sel * 3, 0);
        }

        disp_glob.PasteSumMatrix(&minus_pos_init, 0, 0);

        // paste neighbour not-shared node position vectors
        for (auto neigh_elem_sel = 0; neigh_elem_sel<3; neigh_elem_sel++)
        {
            if (neighbouring_elements[neigh_elem_sel].get() == nullptr)
            {
                disp_glob((neigh_elem_sel + 3) * 3     ) = 0;
                disp_glob((neigh_elem_sel + 3) * 3 + 1 ) = 0;
                disp_glob((neigh_elem_sel + 3) * 3 + 2 ) = 0;
            }
            else
            {
                disp_glob.PasteVector(neighbouring_elements[neigh_elem_sel]->main_nodes[neighbour_node_not_shared[neigh_elem_sel]]->GetPos(), (neigh_elem_sel + 3) * 3, 0);
                disp_glob.PasteSumClippedMatrix(&neighbouring_elements[neigh_elem_sel]->minus_pos_init, neighbour_node_not_shared[neigh_elem_sel] * 3, 0, 3, 1, (neigh_elem_sel + 3) * 3, 0);
            }
                
        }


        

        if (GetNdofs()==18)
        {
            Fi.MatrMultiply(stiffness_matrix, disp_glob);
        }
        else
        {
            Fi_temp.MatrMultiply(stiffness_matrix, disp_glob);

            // paste result into Fi
            Fi.Resize(GetNdofs(), 1);
            Fi.PasteClippedMatrix(&Fi_temp, 0, 0, 9, 1, 0, 0);

            int offset_row = 0;
            for (auto neigh_elem_sel = 0; neigh_elem_sel<3; neigh_elem_sel++)
            {
                if (neighbouring_elements[neigh_elem_sel].get() != nullptr)
                {
                    Fi.PasteClippedMatrix(&Fi_temp, (neigh_elem_sel + 3) * 3, 0, 3, 1, (neigh_elem_sel + 3 - offset_row) * 3, 0);
                }
                else
                    offset_row++;
            }
        }
        

    }

    void SetMaterial(std::shared_ptr<ChMaterialShellTri_3> material)
    {
        m_material = material;
    }

    void SetupInitial(ChSystem* system) override {
        m_material->UpdateConsitutiveMatrices();

        updateBC();

        // Inform the system about which nodes take part to the computation
        std::vector<ChVariables*> vars;
        vars.push_back(&main_nodes[0]->Variables());
        vars.push_back(&main_nodes[1]->Variables());
        vars.push_back(&main_nodes[2]->Variables());

        for (auto neigh_elem_sel = 0; neigh_elem_sel < 3; ++neigh_elem_sel)
        {
            if (neighbouring_elements[neigh_elem_sel].get()!=nullptr)
            {
                vars.push_back(&neighbouring_elements[neigh_elem_sel]->main_nodes[neighbour_node_not_shared[neigh_elem_sel]]->Variables());
            }
        }

        Kmatr.SetVariables(vars);

        // Store the initial global position of nodes;
        minus_pos_init.PasteVector(-main_nodes[0]->GetPos(),0,0);
        minus_pos_init.PasteVector(-main_nodes[1]->GetPos(),0,0);
        minus_pos_init.PasteVector(-main_nodes[2]->GetPos(),0,0);

    }

    /// Gets the number of nodes used by this element.
    int GetNnodes() override { return 3+countNeighbours(); }

    /// Gets the number of coordinates in the field used by the referenced nodes.
    /// This is for example the size (n.of rows/columns) of the local stiffness matrix.
    int GetNdofs() override { return 3 * (3 + countNeighbours()); }

    /// Get the number of coordinates from the n-th node that are used by this element.
    /// Note that this may be different from the value returned by
    ///    GetNodeN(n)->Get_ndof_w();
    int GetNodeNdofs(int n) override { return 3; }

    /// Access the nth node.
    std::shared_ptr<ChNodeFEAbase> GetNodeN(int n) override {
        if (n < 3) // the node belongs to this element
            return main_nodes[n];
        else // the node belongs to a neighbour
        {
            n = n - 3;
            for (auto neigh_elem_sel = 0; neigh_elem_sel<3;neigh_elem_sel++)
            {
                if (neighbouring_elements[neigh_elem_sel]!=nullptr)
                {
                    if (n == 0)
                        return neighbouring_elements[neigh_elem_sel]->main_nodes[neighbour_node_not_shared[neigh_elem_sel]];
                    n--;
                }
            }
            return nullptr;
        }
    }

    void GetStateBlock(ChMatrixDynamic<>& mD) override
    {
        mD.Resize(GetNdofs(), 1);
        mD.PasteVector(this->main_nodes[0]->GetPos(), 0, 0);
        mD.PasteVector(this->main_nodes[1]->GetPos(), 3, 0);
        mD.PasteVector(this->main_nodes[2]->GetPos(), 6, 0);

        int offset_row = 0;
        for (auto neigh_elem_sel =0; neigh_elem_sel<3; neigh_elem_sel++)
        {
            if (neighbouring_elements[neigh_elem_sel] != nullptr)
            {
                mD.PasteVector(this->neighbouring_elements[0]->main_nodes[neighbour_node_not_shared[0]]->GetPos(), 3 * (neigh_elem_sel + 3 - offset_row), 0);
            }
            else
                offset_row++;
        }
    }

    void ComputeMmatrixGlobal(ChMatrix<>& M) override {
        getElementMR(M, 1.0, 0.0);
    }


    //TODO: temporary definition of abstract functions
    void EvaluateSectionDisplacement(const double u,
        const double v,
        const ChMatrix<>& displ,
        ChVector<>& u_displ,
        ChVector<>& u_rotaz) override {};

    void EvaluateSectionFrame(const double u,
        const double v,
        const ChMatrix<>& displ,
        ChVector<>& point,
        ChQuaternion<>& rot) override {};

    void EvaluateSectionPoint(const double u,
        const double v,
        const ChMatrix<>& displ,
        ChVector<>& point) override {};

    void EvaluateSectionVelNorm(double U, double V, ChVector<> &Result) override {};

};

/// @} fea_elements

}  //___end of namespace fea___
}  //___end of namespace chrono___

#endif
