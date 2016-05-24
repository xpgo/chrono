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
// File authors: Dario Mangoni

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
    class ChApiFea ChMaterialShellTri {
    public:
        /// Construct an isotropic material.
        ChMaterialShellTri(double E,    ///< Young's modulus
            double nu,   ///< Poisson ratio
            double rho
        ) :
            m_rho(rho),
            m_E(E),
            m_nu(nu)
        {            
        }

        /// Return the elasticity moduli
        double Get_E() const { return m_E; }
        void Set_E(double E_in) { m_E = E_in; }

        /// Return the Poisson ratio
        double Get_nu() const { return m_nu; }
        void Set_nu(double nu_in) { m_nu = nu_in; }

        /// Return the density
        double Get_rho() const { return m_rho; }
        void Set_rho(double rho_in) { m_rho = rho_in; }

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
            Dbend(2,2) = Dbend_multiplier*(1-m_nu)/2;
        }

        void updateDmembr()
        {
            double Dmembr_multiplier = m_E / (1 - m_nu*m_nu);
            Dbend(0, 0) = Dmembr_multiplier;
            Dbend(0, 1) = m_nu*Dmembr_multiplier;
            Dbend(1, 0) = m_nu*Dmembr_multiplier;
            Dbend(1, 1) = Dmembr_multiplier;
            Dbend(2, 2) = Dmembr_multiplier*(1 - m_nu) / 2;
        }

    private:

        ChMatrixNM<double, 3, 3> Dbend;
        ChMatrixNM<double, 3, 3> Dmembr;
        double m_thickness = 0.1;             ///< thickness
        double m_rho = 7850;                  ///< density
        double m_E = 210e9;                      ///< elasticity moduli
        double m_nu = 0.3;                     ///< Poisson ratio
    };






/// Class for FEA elements of Onate shell type.
/// This element has a linear displacement field.
class ChApiFea ChElementShellTri : public ChElementShell {
protected:
    std::vector<std::shared_ptr<ChNodeFEAxyz>> main_nodes; ///< element nodes

    /// edge versors;
    /// s_loc[1] is of edge 1, that goes from 2 to 3
    /// s_loc[2] is of edge 2, that goes from 3 to 1
    /// s_loc[3] is of edge 3, that goes from 1 to 2
    std::vector<ChMatrixNM<double, 2, 1>> s_loc; 
    std::vector<ChMatrixNM<double, 2, 1>> t_loc; ///< versors tangent to the edge that points inside the main triangle
    ChMatrix33<double> rotGL; ///< rotation matrix from Local to Global frame

private:
    double edge_length[3];
    double element_area = 0;
    double m_thickness = 0;
    ChMatrixNM<double, 18, 18> stiffness_matrix;
    ChMatrixNM<double, 18, 18> mass_matrix;
    ChMatrixNM<double, 18, 18> damping_matrix;
    ChMatrixNM<double, 18, 18> H_matrix;

    

    enum boundary_conditions
    {
        DEFAULT = 0,
        CLAMPED = 1,
        SUPPORTED = 2,
        SYMMETRIC = 3
    } edge_bc[3] = { DEFAULT, DEFAULT, DEFAULT };

    std::vector<std::shared_ptr<ChElementShellTri>> neighbouring_elements; ///< neighbour elements
    /// neighbour_node_not_shared[1] = 3 means that the node of the element attached to the 1st edge has its 3rd node that is not shared with this element
    /// or, that is the same, that the element attached to the 1st edge is attached with its 3rd edge
    int neighbour_node_not_shared[3]; 
    std::shared_ptr<ChMaterialShellTri> m_material;
    double thickness = 0;

    ChMatrix33<double> shape_function;
    ChMatrix33<double> rotLS; ///< rotation matrix from Side to Local frame


    void updateGeometry() {

        // evaluate x' (x local axis) in global coordinates
        ChVector<double> x_loc;
        x_loc = main_nodes[1]->GetPos() - main_nodes[0]->GetPos();
        x_loc.Normalize();
        rotGL.PasteVector(x_loc, 0, 0);
        // evaluate y' (y local axis) in global coordinates
        ChVector<double> y_loc;
        y_loc = main_nodes[2]->GetPos() - main_nodes[0]->GetPos();
        double proj_y_on_x = y_loc.Dot(x_loc);
        y_loc.Sub(y_loc, proj_y_on_x*x_loc );
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
        for (size_t node_sel = 0; node_sel < 3; node_sel++)
        {
            ChVector<double> temp_vect;
            temp_vect = rotLG*main_nodes[node_sel]->GetPos();
            temp_shape(node_sel, 0) = 1;
            temp_shape(node_sel, 1) = temp_vect(1);
            temp_shape(node_sel, 2) = temp_vect(2);
            nodes_pos_loc[node_sel](0, 0) = temp_vect(1);
            nodes_pos_loc[node_sel](1, 0) = temp_vect(2);

        }

        element_area = shape_function.FastInvert(&temp_shape)/2;

        // compute the directions of the edges (not normalized yet)
        s_loc[0] = nodes_pos_loc[2] - nodes_pos_loc[1]; // edge1 goes from 2 to 3
        s_loc[1] = nodes_pos_loc[0] - nodes_pos_loc[2]; // edge2 goes from 3 to 1
        s_loc[2] = nodes_pos_loc[1] - nodes_pos_loc[0]; // edge3 goes from 1 to 2

        // normalize the 's' vectors and build 't' versors
        for (size_t edge_sel = 0; edge_sel < 3; edge_sel++)
        {
            edge_length[edge_sel] = s_loc[edge_sel].NormTwo();
            s_loc[edge_sel].MatrDivScale(edge_length[edge_sel]);
            t_loc[edge_sel](0) = -s_loc[edge_sel](1);
            t_loc[edge_sel](1) = +s_loc[edge_sel](0);
        }


    }


    void updateStructural()
    {
        ///////////////////// Bending stiffness //////////////////////////
        ChMatrixNM<double, 3, 12> Bll; ///< Curvature Matrix (that takes local displacements and returns local curvatures)

        for (size_t edge_sel = 0; edge_sel < 3; edge_sel++)
        {
            ChMatrixNM<double, 2, 12> J;
            
            // theta_s part
            if (!(edge_bc[edge_sel] == CLAMPED || edge_bc[edge_sel] == SYMMETRIC))
            {
                // TODO: do only once at setup step
                // find the edge number respect to the neighbouring element that is shared
                // with the 'edge_sel' of the current element
                if (neighbouring_elements[edge_sel].get()!=nullptr)
                {
                    for (size_t neigh_edge_sel = 0; neigh_edge_sel<3; neigh_edge_sel++)
                    {
                        auto it = std::find(main_nodes.begin(), main_nodes.end(), neighbouring_elements[edge_sel]->main_nodes[neigh_edge_sel]);
                        if (it == main_nodes.end())
                        {
                            neighbour_node_not_shared[edge_sel] = neigh_edge_sel;
                            break;
                        }
                    }

                    // fixed part of theta_s
                    for (size_t col_sel = 0; col_sel<3; col_sel++)
                        J(0, col_sel) = 0.5 * (t_loc[edge_sel](0) * shape_function(1, col_sel) + t_loc[edge_sel](1) * shape_function(2, col_sel));

                    // moving part of theta_s
                    for (size_t col_sel = 0; col_sel<3; col_sel++)
                        J(0, (edge_sel+1)*3+col_sel) = -0.5 * (neighbouring_elements[edge_sel]->t_loc[neighbour_node_not_shared[edge_sel]](0) * neighbouring_elements[edge_sel]->shape_function(1, col_sel)
                                                             + neighbouring_elements[edge_sel]->t_loc[neighbour_node_not_shared[edge_sel]](1) * neighbouring_elements[edge_sel]->shape_function(2, col_sel));
                }
                else
                {
                    // fixed part of theta_s
                    for (size_t col_sel = 0; col_sel<3; col_sel++)
                        J(0, col_sel) = 1 * (t_loc[edge_sel](0) * shape_function(1, col_sel) + t_loc[edge_sel](1) * shape_function(2, col_sel));
                }
            }
                

            // theta_t part
            if ( edge_bc[edge_sel]==CLAMPED || edge_bc[edge_sel] == SUPPORTED )
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
            rot_temp(1, 1) = +s_loc[edge_sel](0)*s_loc[edge_sel](1);
            rot_temp(2, 1) = - rot_temp(0, 0) - rot_temp(1, 0);

            // Update the Curvature Matrix 'Bll'
            ChMatrixNM<double, 3, 12> Bll_temp;
            Bll_temp.MatrMultiply(rot_temp, J);
            Bll_temp.MatrScale(edge_length[edge_sel]);
            Bll.PasteSumMatrix(&Bll_temp,0,0);

        }

        Bll.MatrDivScale(element_area);
        
        // oh man, now we have to compute the rotation matrix for the 'w' displacements; finger crossed...
        ChMatrixNM<double, 18, 12> rotLGw_transp;
        // first: fix the main element rotation part
        rotLGw_transp.PasteClippedMatrix(&rotGL, 0, 2, 3, 1, 0, 0);
        rotLGw_transp.PasteClippedMatrix(&rotGL, 0, 2, 3, 1, 3, 1);
        rotLGw_transp.PasteClippedMatrix(&rotGL, 0, 2, 3, 1, 6, 2);
        for (size_t neigh_elem_sel = 0; neigh_elem_sel < 3; neigh_elem_sel++)
        {
            for (size_t neigh_node_sel = 0; neigh_elem_sel < 3; neigh_elem_sel++)
            {
                // in this loop the algorithm do the following:
                // - takes each node of the selected neighbour element
                // - two nodes of this element are shared with the current element (so are in the main_nodes vector)
                // - in which position are these TWO nodes in the main_nodes vector
                int node_sel;
                for (node_sel = 0; node_sel < 3; node_sel++)
                {
                    if (main_nodes[node_sel]==neighbouring_elements[neigh_elem_sel]->main_nodes[neigh_node_sel])
                    {
                        break;
                    }
                }
                // in node_sel there should be the row of rotLGw_transp in which the rotGL matrix of the current neigh_elem_sel element should be stored
                rotLGw_transp.PasteClippedMatrix(&(neighbouring_elements[neigh_elem_sel]->rotGL), 0, 2, 3, 1, 3*node_sel, 3*(neigh_elem_sel+1)+neigh_node_sel );
            }

            // - the remaining node is NOT shared with the current element
            // - we know which node of the neighbouring element is not shared;
            // - for this we have the 'neighbour_node_not_shared' vector
            rotLGw_transp.PasteClippedMatrix(&neighbouring_elements[neigh_elem_sel]->rotGL, 0, 2, 3, 1, 3 * (neigh_elem_sel+3), 3 * (neigh_elem_sel + 1) + neighbour_node_not_shared[neigh_elem_sel]);
        }

        ChMatrixNM<double, 3, 12> Blg;
        // Update the Curvature Matrix that takes Global displacements and returns Local curvatures
        Blg.MatrMultiplyT(Bll, rotLGw_transp);

        // Update bending stiffness matrix
        ChMatrixNM<double, 18, 18> K_bend;
        K_bend.MatrTMultiply(Blg, m_material->GetConsitutiveMatrixBending());
        K_bend.MatrScale(pow(thickness, 3));
        K_bend.MatrMultiply(K_bend, Blg);

        ///////////////////// Membrane stiffness //////////////////////////
        ChMatrixNM<double, 3, 18> Blg_membr_transp;
        
        ChVector<double> temp = rotGL.ClipVector(0, 0);
        for (size_t col_sel = 0; col_sel < 3; col_sel++)
        {
            Blg_membr_transp.PasteVector(temp*shape_function(1, col_sel), col_sel*3, 0);
            Blg_membr_transp.PasteSumClippedMatrix(&Blg_membr_transp, col_sel * 3, 0, 3, 1, col_sel * 3, 2);
        }

        temp = rotGL.ClipVector(0, 1);
        for (size_t col_sel = 0; col_sel < 3; col_sel++)
        {
            Blg_membr_transp.PasteVector(temp*shape_function(2, col_sel), col_sel * 3, 1);
            Blg_membr_transp.PasteSumClippedMatrix(&Blg_membr_transp, col_sel * 3, 1, 3, 1, col_sel * 3, 2);
        }

        ChMatrixNM<double, 18, 18> K_membr;
        K_membr.MatrMultiply(Blg_membr_transp, m_material->GetConsitutiveMatrixMembrane());
        K_membr.MatrMultiplyT(K_membr, Blg_membr_transp);

        stiffness_matrix.MatrAdd(K_membr, K_bend);
        stiffness_matrix.MatrScale(element_area*thickness);


    }

    void updateBC() //TODO: find other ways to implement boundary conditions
    {
        if (main_nodes[0]->GetFixed() && main_nodes[1]->GetFixed())
            edge_bc[2] == CLAMPED;

        if (main_nodes[1]->GetFixed() && main_nodes[2]->GetFixed())
            edge_bc[0] == CLAMPED;

        if (main_nodes[2]->GetFixed() && main_nodes[0]->GetFixed())
            edge_bc[1] == CLAMPED;
    }

    void updateK()
    {
        updateGeometry();
        updateStructural();
    }

    void updateM()
    {
        mass = element_area * m_thickness * m_material->Get_rho();
        for (size_t diag_sel = 0; diag_sel < 3; diag_sel++)
        {
            mass_matrix(diag_sel * 3    , diag_sel * 3    ) = mass / 3;
            mass_matrix(diag_sel * 3 + 1, diag_sel * 3 + 1) = mass / 3;
            mass_matrix(diag_sel * 3 + 2, diag_sel * 3 + 2) = mass / 3;
        }
    }

    void updateR()
    {
        double damp_temp = 1e-3;
        for (size_t diag_sel = 0; diag_sel < 18; diag_sel++)
        {
            damping_matrix(diag_sel, diag_sel) = damp_temp;
        }
    }




public:
    ChElementShellTri() {
        main_nodes.resize(3);
        neighbouring_elements.resize(3);
        s_loc.resize(3);
        t_loc.resize(3);
    }

    /// Return the thickness
    double Get_thickness() const { return m_thickness; }
    void Set_thickness(double thickness_in) { m_thickness = thickness_in; }

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
         
        for (size_t elem_sel = 0; elem_sel < mesh->GetNelements(); elem_sel++)
        {
            int common_nodes[3] = {-1,-1,-1};
            int common_nodes_count = 0;
            for (size_t node_sel = 0; node_sel < 3; node_sel++)
            {
                for (size_t main_node_sel = 0; main_node_sel < 3; main_node_sel++)
                {
                    if (mesh->GetElement(elem_sel)->GetNodeN(node_sel) == main_nodes[main_node_sel])
                    {
                        // the element has at least a node in common

                        if (mesh->GetElement(elem_sel).get() == this) // the element is being compared with itself!
                            break;

                        common_nodes[main_node_sel] = node_sel;
                        common_nodes_count++;

                    }
                }
            }
            if (common_nodes_count > 1)
            {
                // find the non shared node
                for (size_t main_node_sel = 0; main_node_sel < 3; main_node_sel++)
                {
                    if (common_nodes[main_node_sel] == -1) // look for not shared node within the main_nodes
                    {
                        // 'main_node_sel' points at the node that it is not shared;
                        // 'main_node_sel' is also the edge on which the element is attached
                        neighbouring_elements[main_node_sel] = std::dynamic_pointer_cast<ChElementShellTri>(mesh->GetElement(elem_sel));
                        // look which node of the neighbour is not linked to the current element
                        // if it is not linked the number won't be in common_nodes[...]
                        for (size_t missing_number = 0; missing_number<3; ++missing_number)  // look for node number 0, 1 and 2
                        {
                            size_t main_node_sel2;
                            for (main_node_sel2=0; main_node_sel2<3; ++main_node_sel2)
                            {
                                // look if in common_nodes[main_node_sel2] there is the 'missing_number-th' node of the selected element
                                if (common_nodes[main_node_sel2] == missing_number)
                                    break; // if there is then 'missing_number' is not missing;
                            }

                            // if 'missing_number' is not missing then main_node_sel2 will be less than 3
                            if (main_node_sel2 == 3) 
                            {
                                neighbour_node_not_shared[main_node_sel] = missing_number;
                                break;
                            }
                                
                        }
                        
                    }
                        
                }
            }

        }
    }
    

    void SetBC(int BC_on_edge1, int BC_on_edge2, int BC_on_edge3)
    {
        edge_bc[0] = static_cast<boundary_conditions>(BC_on_edge1);
        edge_bc[1] = static_cast<boundary_conditions>(BC_on_edge2);
        edge_bc[3] = static_cast<boundary_conditions>(BC_on_edge3);
    }

    ChMatrixNM<double, 3, 3>& GetShapeFunction() { return shape_function; }

    void ComputeKRMmatricesGlobal(ChMatrix<>& H, double Kfactor, double Rfactor, double Mfactor) override {
        H_matrix = stiffness_matrix;
        H_matrix.MatrScale(Kfactor);
        
        ChMatrixNM<double, 18, 18> temp_mat;
        temp_mat = mass_matrix;
        temp_mat.MatrScale(Mfactor);
        H_matrix.PasteSumMatrix(&temp_mat, 0, 0);

        temp_mat = damping_matrix;
        temp_mat.MatrScale(Rfactor);
        H_matrix.PasteSumMatrix(&temp_mat, 0, 0);
        
    }

    void ComputeInternalForces(ChMatrixDynamic<>& Fi) override{
        Fi.Resize(18,1); // also clears the vector
        ChMatrixNM<double, 18, 1> noder_pos_glob;
        ChMatrixNM<double, 3, 1> temp;
        for (size_t node_sel = 0; node_sel<3; node_sel++)
        {
            noder_pos_glob.PasteVector(main_nodes[node_sel]->GetPos(), node_sel * 3, 0);
        }

        Fi.MatrMultiply(stiffness_matrix, noder_pos_glob);

    }

    void SetMaterial(std::shared_ptr<ChMaterialShellTri> material)
    {
        m_material = material;
    }

    void SetupInitial(ChSystem* system) override {
        m_material->UpdateConsitutiveMatrices();
    }

    /// Gets the number of nodes used by this element.
    int GetNnodes() override { return 3; }

    /// Gets the number of coordinates in the field used by the referenced nodes.
    /// This is for example the size (n.of rows/columns) of the local stiffness matrix.
    int GetNdofs() override { return 18; }

    /// Get the number of coordinates from the n-th node that are used by this element.
    /// Note that this may be different from the value returned by
    ///    GetNodeN(n)->Get_ndof_w();
    int GetNodeNdofs(int n) override { return 3; }

    /// Access the nth node.
    std::shared_ptr<ChNodeFEAbase> GetNodeN(int n) override { return main_nodes[n]; }

    void GetStateBlock(ChMatrixDynamic<>& mD) override
    {
        mD.Resize(18, 1);
        //TODO: what does this function?
    };

    void ComputeMmatrixGlobal(ChMatrix<>& M) override {
        M = static_cast<ChMatrix<double>>(mass_matrix);
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
