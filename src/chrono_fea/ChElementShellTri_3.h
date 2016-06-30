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
#include <array>

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
public: // TODO: IT IS PUBLIC ONLY FOR DEBUG PURPOSE!!!!
    std::array<std::shared_ptr<ChNodeFEAxyz>,6> all_nodes; ///< it is just [main_nodes, neigh_nodes]

    static constexpr std::array<std::array<int, 3>, 4> numeration = { {{0,1,2},{3,2,1},{4,0,2},{5,1,0}} };
    static constexpr std::array<std::array<int, 4>, 4> edge_num = { {{1,2,0,1},{2,1,3,2},{0,2,4,0},{1,0,5,1}} };

    double element_area0 = -1;
    double thickness = -1;
    ChMatrixNM<double, 18, 18> stiffness_matrix;

    enum boundary_conditions
    {
        DEFAULT = 0,
        CLAMPED = 1,
        SUPPORTED = 2,
        SYMMETRIC = 3
    } edge_bc[3] = { DEFAULT, DEFAULT, DEFAULT };

    std::array<std::shared_ptr<ChElementShellTri_3>,3> neighbouring_elements; ///< neighbour elements

    std::shared_ptr<ChMaterialShellTri_3> m_material;


    double GetArea0()
    {
        if (element_area0 < 0)
            initializeElement();

        return element_area0;
    }

    //inline void GetEdgeVector0(ChVector<double>& vect, int edge_sel, int elem_sel, ChVector<double> (ChNodeFEAxyz::*fun)())
    //{
    //    vect.Sub(all_nodes[GetNodeOfEdge(edge_sel, 1, elem_sel)]->(*fun)(), all_nodes[GetNodeOfEdge(edge_sel, 0, elem_sel)]->(*fun)());
    //}

    inline void GetEdgeVector0(ChVector<double>& vect, int edge_sel, int elem_sel = 0) const
    {
        vect.Sub(all_nodes[GetNodeOfEdge(edge_sel, 1, elem_sel)]->GetX0(), all_nodes[GetNodeOfEdge(edge_sel, 0, elem_sel)]->GetX0());
    }

    inline void GetEdgeVector(ChVector<double>& vect, int edge_sel, int elem_sel = 0) const
    {
        vect.Sub(all_nodes[GetNodeOfEdge(edge_sel, 1, elem_sel)]->GetPos(), all_nodes[GetNodeOfEdge(edge_sel, 0, elem_sel)]->GetPos());
    }

    ChVector<double> GetEdgeVector(int edge_sel, int elem_sel = 0) const
    {
        ChVector<double> vect;
        vect.Sub(all_nodes[GetNodeOfEdge(edge_sel, 1, elem_sel)]->GetPos(), all_nodes[GetNodeOfEdge(edge_sel, 0, elem_sel)]->GetPos());
        return vect;
    }
    
    inline double GetEdgeVersor0(ChVector<double>& vers, int edge_sel, int elem_sel = 0) const
    {
        GetEdgeVector0(vers, edge_sel, elem_sel);
        double norm = vers.Length();
        vers.Scale(1/norm);
        return norm;
    }

    static int inline GetNodeOfEdge(int edge_sel, int node_sel, int elem_sel = 0)
    {
        return edge_num[elem_sel][edge_sel + node_sel];
    }

    double GetHeight(int edge_sel, int elem_sel) const
    {
        ChVector<double> v_edge_sel;
        ChVector<double> v_edge_diff;
        GetEdgeVector(v_edge_sel, edge_sel, elem_sel);
        GetEdgeVector(v_edge_diff, GetNodeOfEdge(edge_sel,0,elem_sel), elem_sel);
        v_edge_sel.Scale(v_edge_sel.Dot(v_edge_diff));
        v_edge_sel+=v_edge_diff;
        return v_edge_sel.Length();
    }

    std::array<ChVector<double>, 3> edge_versors0;
    std::array<double, 3> edge_length0;
    std::array<ChVector<double>, 3> ip_normal_mod; // in-plane normals by a factor of edge_len/2/A
    ChMatrixNM<double, 3, 6> J_source;

    void initializeElement()
    {
        if (element_area0 > 0)
            return;

        /* Stores variables to be used to compute the curvature matrix */
        //store edge versors
        for (auto edge_sel = 0; edge_sel<3; ++edge_sel)
        {
            edge_length0[edge_sel] = GetEdgeVersor0(edge_versors0[edge_sel], edge_sel, 0);
        }

        //compute area of main element
        ChVector<double> vect_temp;    vect_temp.Cross(edge_versors0[0], edge_versors0[1]);
        element_area0 = edge_length0[0]* edge_length0[1]*vect_temp.Length()/2; //TODO: check if the area equation is correct


        //normal in-plane versors of the main element M (glob ref)
        //project the triangle on the XY plane
        //compute edge versors, rotate -90° around Z (outward normal)
        for (auto edge_sel = 0; edge_sel<3; ++edge_sel)
        {
            double norm_temp = edge_versors0[edge_sel](0)*edge_versors0[edge_sel](0) + edge_versors0[edge_sel](1)*edge_versors0[edge_sel](1);
            ip_normal_mod[edge_sel](0) =  edge_versors0[edge_sel](1) / norm_temp * edge_length0[edge_sel] / 2 / element_area0;
            ip_normal_mod[edge_sel](1) = -edge_versors0[edge_sel](0) / norm_temp * edge_length0[edge_sel] / 2 / element_area0;
        }

        // compute scale factors
        for (auto edge_sel = 0; edge_sel<3; ++edge_sel)
        {
            for (auto col_sel = 0; col_sel<3; ++col_sel)
            {
                J_source(edge_sel, col_sel) = edge_length0[col_sel] * ChVector<double>::Dot(edge_versors0[col_sel], edge_versors0[edge_sel]) / 2 / element_area0;
                if (neighbouring_elements[edge_sel])
                    J_source(edge_sel, col_sel + 3) = -ChVector<double>::Dot(GetEdgeVector(col_sel, edge_sel), edge_versors0[edge_sel]) / 2 / neighbouring_elements[edge_sel]->GetArea0();
            }
        }
        /* ********************************************************** */

    }



    int countNeighbours() const
    {
        int count = 0;
        for (auto node_sel : all_nodes)
        {
            if (node_sel)
            {
                count++;
            }
        }
        return count;
    }

    void updateBC() //TODO: find other ways to implement boundary conditions
    {
        if (all_nodes[0]->GetFixed() && all_nodes[1]->GetFixed())
            edge_bc[2] = CLAMPED;

        if (all_nodes[1]->GetFixed() && all_nodes[2]->GetFixed())
            edge_bc[0] = CLAMPED;

        if (all_nodes[2]->GetFixed() && all_nodes[0]->GetFixed())
            edge_bc[1] = CLAMPED;
    }

    void updateElementMass(int node_sel)
    {
        mass = element_area0 * thickness * m_material->Get_rho();
        for (auto diag_sel = 0; diag_sel < 3; diag_sel++)
        {
            all_nodes[diag_sel]->SetMass(mass/3);
        }
    }

    void getNodeM(ChMatrix<double>& node_mass, int node_sel, double factor = 1.0)
    {
        if (node_mass.GetColumns() != 3 || node_mass.GetRows() != 3)
            node_mass.Resize(3, 3);

        updateElementMass(node_sel);
        node_mass.FillDiag(all_nodes[node_sel]->GetMass()*factor);

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


public:
    ChElementShellTri_3()
    {
    }

    /// Return the thickness
    double Get_thickness() const { return thickness; }
    void Set_thickness(double thickness_in) { thickness = thickness_in; }

    /// Sets only the nodes of the main element; the neighbouring will be added by SetupInitial()
    void SetNodes(std::shared_ptr<ChNodeFEAxyz> node0,
                  std::shared_ptr<ChNodeFEAxyz> node1,
                  std::shared_ptr<ChNodeFEAxyz> node2)
    {
        all_nodes[0] = node0;
        all_nodes[1] = node1;
        all_nodes[2] = node2;
    }

    /// Inform the current element of which nodes are in its neighbour
    void SetNeighbouringNodes(std::shared_ptr<ChNodeFEAxyz> node3,
        std::shared_ptr<ChNodeFEAxyz> node4,
        std::shared_ptr<ChNodeFEAxyz> node5)
    {
        all_nodes[3] = node3;
        all_nodes[4] = node4;
        all_nodes[5] = node5;
    }

    
    void Update() override {
        /* ************** Compute the curvature matrix ************** */
        // compute r_i^M
        std::array<double, 3> stiff_ratio;
        for (auto edge_sel = 0; edge_sel<3; ++edge_sel)
        {
            if (neighbouring_elements[edge_sel])
                stiff_ratio[edge_sel] = 1 / (1 + (m_material->Get_E() *pow(thickness, 3) * GetHeight(0, edge_sel)) / (neighbouring_elements[edge_sel]->m_material->Get_E() *pow(neighbouring_elements[edge_sel]->thickness, 3) * GetHeight(edge_sel, 0)));
            else
                if (edge_bc[edge_sel] == CLAMPED) // only clamped edge BC supported
                    stiff_ratio[edge_sel] = 1;
                else
                    stiff_ratio[edge_sel] = 0;

        }

        // compute normal to plane
        std::array<ChVector<double>, 4> t;
        for (auto elem_sel = 0; elem_sel<4; ++elem_sel)
        {
            if (all_nodes[elem_sel+2]) // check if neighbour exists before trying to compute area
            {
                t[elem_sel].Cross(GetEdgeVector(0, elem_sel), GetEdgeVector(1, elem_sel));
                t[elem_sel].Normalize();
            }
            
        }

        // compute J
        ChMatrixNM<double, 3, 18> B_bend;
        for (auto edge_sel = 0; edge_sel<3; ++edge_sel)
        {
            for (auto col_sel = 0; col_sel<3; ++col_sel)
            {
                B_bend.PasteSumVector(ip_normal_mod[edge_sel] * (stiff_ratio[edge_sel] * t[0](col_sel)* J_source(edge_sel, col_sel)), 0, col_sel * 3);
                B_bend.PasteSumVector(ip_normal_mod[edge_sel] * (stiff_ratio[edge_sel] * t[0](col_sel)* J_source(edge_sel, col_sel)), 0, col_sel * 3 + 1);
                B_bend.PasteSumVector(ip_normal_mod[edge_sel] * (stiff_ratio[edge_sel] * t[0](col_sel)* J_source(edge_sel, col_sel)), 0, col_sel * 3 + 2);

                if (all_nodes[edge_sel + 3]) // check if neighbour exists before trying to compute area
                {
                    B_bend.PasteSumVector(ip_normal_mod[edge_sel] * (stiff_ratio[edge_sel] * t[edge_sel + 1](col_sel)* J_source(edge_sel, col_sel)), 0, col_sel * 3);
                    B_bend.PasteSumVector(ip_normal_mod[edge_sel] * (stiff_ratio[edge_sel] * t[edge_sel + 1](col_sel)* J_source(edge_sel, col_sel)), 0, col_sel * 3 + 1);
                    B_bend.PasteSumVector(ip_normal_mod[edge_sel] * (stiff_ratio[edge_sel] * t[edge_sel + 1](col_sel)* J_source(edge_sel, col_sel)), 0, col_sel * 3 + 2);
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


    /// Finds the neighbouring elements and nodes
    /// To be called after that ALL the elements in the mesh have their main nodes set.
    void UpdateConnectivity(std::shared_ptr<ChMesh> mesh)
    {
        for (unsigned int elem_sel = 0; elem_sel<mesh->GetNelements(); ++elem_sel) // pick an element in the list that is candidate to be a neighbour
        {
            int shared_nodes_counter = 0;
            ChNodeFEAbase* not_shared_node_ptr = nullptr;
            for (auto neigh_node_sel = 0; neigh_node_sel<3; ++neigh_node_sel) // pick a node of the candidate element
            {
                int main_node_sel;
                for (main_node_sel = 0; main_node_sel<3; ++main_node_sel) // pick a node of the current element
                {
                    if (mesh->GetElement(elem_sel)->GetNodeN(neigh_node_sel) == all_nodes[main_node_sel]) // if the two nodes picked are different...
                    {
                        shared_nodes_counter++; // counts how many shared nodes has been found
                        break;
                    }
                }

                if (main_node_sel==3)
                {
                    if (!not_shared_node_ptr)
                        not_shared_node_ptr = mesh->GetElement(elem_sel)->GetNodeN(neigh_node_sel).get();
                    else
                        break;
                }
            }

            if (!not_shared_node_ptr && shared_nodes_counter==2) // exclude the case that the candidate and current elements coincide;
            {

                // on which edge the just-found neighbour is attached?
                int not_shared_node_sel; // the not_shared_node_sel will point to the node in all_nodes that is not shared
                // --> it will be also the number of the edge on which it is attached
                for (not_shared_node_sel = 0; not_shared_node_sel<3; ++not_shared_node_sel) 
                {
                    if (all_nodes[not_shared_node_sel].get() != not_shared_node_ptr)
                    {
                        break;
                    }
                }

                all_nodes[3 + not_shared_node_sel] = std::make_shared<ChNodeFEAxyz>(*dynamic_cast<ChNodeFEAxyz*>(not_shared_node_ptr));
                neighbouring_elements[not_shared_node_sel] = std::dynamic_pointer_cast<ChElementShellTri_3>(mesh->GetElement(elem_sel));
            }

        }

    }


    void SetBC(int BC_on_edge1, int BC_on_edge2, int BC_on_edge3)
    {
        edge_bc[0] = static_cast<boundary_conditions>(BC_on_edge1);
        edge_bc[1] = static_cast<boundary_conditions>(BC_on_edge2);
        edge_bc[2] = static_cast<boundary_conditions>(BC_on_edge3);
    }



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

        // TODO: it shouldn't be necessary since K it is updated during Update() call
        Update();


        // paste main_nodes position vectors
        for (auto node_sel = 0; node_sel<3; node_sel++){
            disp_glob.PasteVector( all_nodes[node_sel]->GetPos() - all_nodes[node_sel]->GetX0(), node_sel * 3, 0);
        }

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
                //disp_glob.PasteVector(neighbouring_elements[neigh_elem_sel]->all_nodes[neighbour_node_not_shared[neigh_elem_sel]]->GetPos()
                //                      - neighbouring_elements[neigh_elem_sel]->all_nodes[neighbour_node_not_shared[neigh_elem_sel]]->GetX0(), (neigh_elem_sel + 3) * 3, 0);
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

    /// Gets the number of nodes used by this element.
    int GetNnodes() override { return 3+countNeighbours(); }

    /// Gets the number of coordinates in the field used by the referenced nodes.
    /// This is for example the size (n.of rows/columns) of the local stiffness matrix.
    int GetNdofs() override { return 3 * (3 + countNeighbours()); }

    /// Get the number of coordinates from the n-th node that are used by this element.
    /// Note that this may be different from the value returned by
    ///    GetNodeN(n)->Get_ndof_w();
    int GetNodeNdofs(int n) override { return 3; }

    /// Access the nth node (it retrieves also the nodes that actually belongs to neighbours)
    std::shared_ptr<ChNodeFEAbase> GetNodeN(int n) override {
        return all_nodes[n];
    }

    void GetStateBlock(ChMatrixDynamic<>& mD) override
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

    void ComputeMmatrixGlobal(ChMatrix<>& M) override {
        getElementMR(M, 1.0, 0.0);
    }


    //TODO: temporary definition of abstract functions. Needed since we are inheriting from ChElementShell
    void EvaluateSectionDisplacement(const double u,
        const double v,
        const ChMatrix<>& displ,
        ChVector<>& u_displ,
        ChVector<>& u_rotaz) override {}

    void EvaluateSectionFrame(const double u,
        const double v,
        const ChMatrix<>& displ,
        ChVector<>& point,
        ChQuaternion<>& rot) override {}

    void EvaluateSectionPoint(const double u,
        const double v,
        const ChMatrix<>& displ,
        ChVector<>& point) override {}

    void EvaluateSectionVelNorm(double U, double V, ChVector<> &Result) override {}

};

/// @} fea_elements

}  //___end of namespace fea___
}  //___end of namespace chrono___

#endif
