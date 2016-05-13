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

namespace chrono {
namespace fea {

/// @addtogroup fea_elements
/// @{

/// Material definition.
/// This class implements material properties for a layer from the Reissner theory,
/// see Morandini, Masarati "Implementation and Validation of a 4-node shell element"
    class ChApiFea ChMaterialShellTri {
    public:
        /// Construct an isotropic material.
        ChMaterialShellTri(double thickness, ///< thickness
            double E,    ///< Young's modulus
            double nu,   ///< Poisson ratio
        );

        /// Return the thickness
        double Get_thickness() const { return m_thickness; }

        /// Return the elasticity moduli
        double Get_E() const { return m_E; }
        /// Return the Poisson ratio
        double Get_nu() const { return m_nu; }

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
        double m_thickness;             ///< thickness
        double m_rho;                  ///< density
        double m_E;                      ///< elasticity moduli
        double m_nu;                     ///< Poisson ratio
    };

/// Class for FEA elements of Onate shell type.
/// This element has a linear displacement field.
class ChApiFea ChElementShellTri {
private:
    friend class ChPatchShellTri;
    std::vector<std::shared_ptr<ChNodeFEAxyz>> m_nodes; ///< element nodes
    std::shared_ptr<ChMaterialShellTri> m_material;

    ChMatrix33<double> shape_function_coefficients;
    
    ChMatrix33<double> rotLS; ///< rotation matrix from Side to Local frame

    void GetElementData() {
        // Rotation matrix from Local to Global frame: rotGL
        ChMatrix33<double> rotGL;
        // evaluate x' (x local axis) in global coordinates
        ChVector<double> x_loc;
        x_loc = m_nodes[1]->GetPos() - m_nodes[0]->GetPos();
        x_loc.Normalize();
        rotGL.PasteVector(x_loc, 0, 0);
        // evaluate y' (y local axis) in global coordinates
        ChVector<double> y_loc;
        y_loc = m_nodes[2]->GetPos() - m_nodes[0]->GetPos();
        double proj_y_on_x = y_loc.Dot(x_loc);
        y_loc.Sub(y_loc, proj_y_on_x*x_loc );
        y_loc.Normalize();
        rotGL.PasteVector(y_loc, 0, 1);
        // evaluate z' (z local axis) in global coordinates
        ChVector<double> z_loc;
        z_loc.Cross(x_loc, y_loc);
        rotGL.PasteVector(z_loc, 0, 2);

        ChMatrix33<double> rotLG;
        double rotLG_det = rotGL.FastInvert(&rotLG);
        // Local position of nodes
        std::vector<std::shared_ptr<ChVector<double>>> nodes_loc;
        nodes_loc[0] = rotLG.

        //ChMatrix33<double> temp;
        //ChVector<double> pos_temp;
        //temp(0, 0) = 1; temp(1, 0) = 1; temp(2, 0) = 1;
        //pos_temp = m_nodes[0]->GetPos();
        //temp(0, 1) = pos_temp()
    }

    void SetNeighbours() { // when could I call it?

    }

public:
    ChElementShellTri() {

    }

    void SetNodes(std::shared_ptr<ChNodeFEAxyz> nodeA,
                  std::shared_ptr<ChNodeFEAxyz> nodeB,
                  std::shared_ptr<ChNodeFEAxyz> nodeC)
    {
        m_nodes[0] = nodeA;
        m_nodes[1] = nodeB;
        m_nodes[2] = nodeC;
    };

    /// Fills the N shape function matrix.
    void ShapeFunctions(ChMatrix<>& N, double x, double y, double z) {

    };

};

/// @} fea_elements

}  //___end of namespace fea___
}  //___end of namespace chrono___

#endif
