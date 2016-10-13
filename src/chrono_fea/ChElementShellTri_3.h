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
            m_density(rho),
            m_YoungModulus(E),
            m_PoissonRatio(nu)
        {
            UpdateConsitutiveMatrices();
        }

        /// Return the elasticity moduli
        double Get_E() const { return m_YoungModulus; }
        void Set_E(double E_in) { m_YoungModulus = E_in; UpdateConsitutiveMatrices(); }

        /// Return the Poisson ratio
        double Get_nu() const { return m_PoissonRatio; }
        void Set_nu(double nu_in) { m_PoissonRatio = nu_in; UpdateConsitutiveMatrices();
        }

        /// Return the density
        double Get_rho() const { return m_density; }
        void Set_rho(double rho_in) { m_density = rho_in; UpdateConsitutiveMatrices();
        }

        void UpdateConsitutiveMatrices()
        {
            constitutive_matrix(0, 0) = 1.0f;
            constitutive_matrix(0, 1) = m_PoissonRatio;
            constitutive_matrix(1, 0) = m_PoissonRatio;
            constitutive_matrix(1, 1) = 1.0f;
            constitutive_matrix(2, 2) = (1.0f - m_PoissonRatio) / 2.0f;
            constitutive_matrix *= m_YoungModulus / (1.0f - m_PoissonRatio*m_PoissonRatio);
        }

        const ChMatrixNM<double, 3, 3>& GetConsitutiveMatrix() const { return constitutive_matrix; };

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

    private:

        ChMatrixNM<double, 3, 3> constitutive_matrix;

        double m_density = 7850;    ///< density
        double m_YoungModulus = 210e9; ///< elasticity moduli
        double m_PoissonRatio = 0.3;  ///< Poisson ratio
    };





/// Class for FEA elements of Onate shell type.
/// This element has a linear displacement field,
/// so the curvature is computed using position of adjacent out-of-element nodes.
class ChApiFea ChElementShellTri_3 : public ChElementShell {
private:
public: // TODO: PUBLIC ONLY FOR DEBUG PURPOSE!!!!

    std::array<std::shared_ptr<ChNodeFEAxyz>,6> all_nodes; ///< the first three elements are the nodes of the current element [main_nodes]; last three are the neighbours' node [neigh_nodes]
    std::array<std::shared_ptr<ChElementShellTri_3>, 3> neighbouring_elements; ///< neighbour elements

    // Connectivity variables
    /// each element owns three 'main_nodes' that are, at the same time, 'neigh_nodes' for other elements. Now:
    /// - each row ('i') refers to an element; (0: main, 1: element0, 2: element1, 3: element2)
    /// - the j-th column refers to the node number in the neighbours' local numeration
    /// - neigh_to_local_numbering[i][j] tells in which position in 'all_nodes' we can find the j-th node of the i-th element
    static constexpr std::array<std::array<int, 3>, 4> neigh_to_local_numbering = { {{0,1,2},{2,1,3},{0,2,4},{1,0,5}} };

    /// each row ('i') refers to an element; (0: main, 1: element0, 2: element1, 3: element2)
    /// the n-th edge goes of element 'i' goes from edge_num[i][n] to edge_num[i][n+1]
    static constexpr std::array<std::array<int, 4>, 4> edge_num = { {{1,2,0,1},{2,1,3,2},{0,2,4,0},{1,0,5,1}} };

    // Geometric variables
    double thickness = -1;
    double thickness0 = 0.025;

    std::shared_ptr<ChMaterialShellTri_3> m_material;

    // internal use
    ChMatrixNM<double, 18, 18> stiffness_matrix;

    enum class boundary_conditions
    {
        DEFAULT,
        CLAMPED,
        SUPPORTED,
        SYMMETRIC
    } edge_bc[3] = {boundary_conditions::DEFAULT, boundary_conditions::DEFAULT, boundary_conditions::DEFAULT };


    //TODO: many of these variables below are only for debug purpose
    

    ChVector<double> z_versor0;
    std::array<double, 3> edge_length0 = { -1,-1,-1 };
    ChMatrix33<double> rotGL0;
    ChMatrixNM<double, 2, 3> edge_normal_vers0;
    ChMatrixNM<double, 3, 2> gradient_shape_function0;
    ChMatrixNM<double, 3, 2> gradient_local0;
    ChMatrixNM<double, 3, 1> gradient_side_n0;
    ChMatrix33<double> c_proj0;
    double element_area0 = -1;

    std::array<ChMatrixNM<double,3,1>,3> L_block0;


    size_t update_counter = 0;


    static int getLocalNodeNumber(int element, int neighbour_numbering);
    // internal functions
    int countNeighbours() const; ///< counts how many neighbours are actually connected to this element
    void updateBC(); //TODO: find other ways to implement boundary conditions
    void updateElementMass(int node_sel);
    void getNodeM(ChMatrix<double>& node_mass, int node_sel, double factor = 1.0);
    void getNodeK(ChMatrix<double>& node_damp, int node_sel, double factor = 1.0);
    void getElementMR(ChMatrix<>& H, double Mfactor, double Kfactor);
    void getElementData(const ChVector<>& P1_in, const ChVector<>& P2_in, const ChVector<>& P3_in,
                        std::array<double, 3>* edge_length,
                        ChVector<>* z_versor,
                        double* area,
                        ChMatrix33<>* rotGL,
                        ChMatrix<>* gradient_shape_function,
                        ChMatrix<>* gradient_local,
                        ChMatrix<>* gradient_side_n,
                        ChMatrix<>* edge_normal_vers,
                        ChMatrix<>* c_proj);

public:
    ChElementShellTri_3(){}
    virtual ~ChElementShellTri_3(){}

    /// Return the thickness
    double Get_thickness() const { return thickness; }
    void Set_thickness(double thickness_in) { thickness = thickness_in; }

    void Update() override;

    /// Finds the neighbouring elements and nodes
    /// To be called after that ALL the elements in the mesh have their main nodes set.
    void UpdateConnectivity(std::shared_ptr<ChMesh> mesh);
    void SetBC(int BC_on_edge1, int BC_on_edge2, int BC_on_edge3);

    void ComputeKRMmatricesGlobal(ChMatrix<>& H, double Kfactor, double Rfactor, double Mfactor) override;
    void ComputeInternalForces(ChMatrixDynamic<>& Fi) override;
    void SetMaterial(std::shared_ptr<ChMaterialShellTri_3> material) { m_material = material; }

    void SetupInitial(ChSystem* system) override;

    /// Gets the number of nodes used by this element.
    int GetNnodes() override { return 3 + countNeighbours(); }

    /// Gets the number of coordinates in the field used by the referenced nodes.
    /// This is for example the size (n.of rows/columns) of the local stiffness matrix.
    int GetNdofs() override { return 3 * (3 + countNeighbours()); }

    /// Get the number of coordinates from the n-th node that are used by this element.
    /// Note that this may be different from the value returned by
    ///    GetNodeN(n)->Get_ndof_w();
    int GetNodeNdofs(int n) override { return 3; }

    /// Access the nth node (it retrieves also the nodes that actually belongs to neighbours)
    std::shared_ptr<ChNodeFEAbase> GetNodeN(int n) override;

    void GetStateBlock(ChMatrixDynamic<>& mD) override;

    void ComputeMmatrixGlobal(ChMatrix<>& M) override;

    /// Sets only the nodes of the main element; the neighbouring will be added by SetupInitial()
    void SetNodes(std::shared_ptr<ChNodeFEAxyz> node0,
        std::shared_ptr<ChNodeFEAxyz> node1,
        std::shared_ptr<ChNodeFEAxyz> node2);

    /// Sets only the nodes of the main element; the neighbouring will be added by SetupInitial()
    void SetNodes(std::shared_ptr<ChNodeFEAxyz> node0,
        std::shared_ptr<ChNodeFEAxyz> node1,
        std::shared_ptr<ChNodeFEAxyz> node2,
        std::shared_ptr<ChNodeFEAxyz> node3,
        std::shared_ptr<ChNodeFEAxyz> node4,
        std::shared_ptr<ChNodeFEAxyz> node5);

    /// Inform the current element of which nodes are in its neighbourhood
    void SetNeighbouringNodes(std::shared_ptr<ChNodeFEAxyz> node3,
        std::shared_ptr<ChNodeFEAxyz> node4,
        std::shared_ptr<ChNodeFEAxyz> node5);


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
