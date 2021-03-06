// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Rainer Gericke
// =============================================================================
//
// HMMWV TMeasy tire subsystem
//
// =============================================================================

#ifndef HMMWV_TMEASY_TIRE_H
#define HMMWV_TMEASY_TIRE_H

#include "chrono/assets/ChTriangleMeshShape.h"

#include "chrono_vehicle/wheeled_vehicle/tire/ChTMeasyTire.h"

#include "chrono_models/ChApiModels.h"

namespace chrono {
namespace vehicle {
namespace hmmwv {

/// @addtogroup vehicle_models_wvp
/// @{

/// TMeasy tire model for the WVP vehicle.
class CH_MODELS_API HMMWV_TMeasyTire : public ChTMeasyTire {
  public:
    HMMWV_TMeasyTire(const std::string& name);
    ~HMMWV_TMeasyTire() {}

    virtual double GetVisualizationWidth() const override { return m_width; }

    virtual void SetTMeasyParams() override;
    virtual double GetMass() const override { return m_mass; }
    virtual ChVector<> GetInertia() const override { return m_inertia; }

    virtual void AddVisualizationAssets(VisualizationType vis) override;
    virtual void RemoveVisualizationAssets() override final;

    void GenerateCharacteristicPlots(const std::string& dirname);

  private:
    static const std::string m_meshName;
    static const std::string m_meshFile;
    static const double m_mass;
    static const ChVector<> m_inertia;

    std::shared_ptr<ChTriangleMeshShape> m_trimesh_shape;
    ChFunction_Recorder m_stiffnessMap;
};

/// @} vehicle_models_hmmwv

}  // end namespace hmmwv
}  // end namespace vehicle
}  // end namespace chrono

#endif
