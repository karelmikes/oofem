/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef misesmattension_h
#define misesmattension_h

#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/structuralms.h"
#include "sm/Materials/isolinearelasticmaterial.h"
#include "dictionary.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "scalarfunction.h"

//#include "linearelasticmaterial.h"
//#include "floatmatrixf.h"
//#include "floatarrayf.h"
#include "sm/Materials/qcmaterialextensioninterface.h"

///@name Input fields for MisesMatTension
//@{
#define _IFT_MisesMatTension_Name "MisesMatTension"
#define _IFT_MisesMatTension_sig0 "sig0"
#define _IFT_MisesMatTension_h "h"
#define _IFT_MisesMatTension_omega_crit "omega_crit"
#define _IFT_MisesMatTension_a "a"
#define _IFT_MisesMatTension_yieldTol "yieldtol"
//@}

namespace oofem {
class GaussPoint;

/**
 * This class implements an isotropic linear elastic material in a finite element problem.
 * For large deformation analysis it becomes the St. Venant-Kirchoff hyperelasticity model.
 *
 * Tasks:
 * - Returning standard material stiffness matrix for 3d-case.
 *   according to current state determined by using data stored
 *   in Gausspoint.
 * - methods give2dPlaneStressMtrx, givePlaneStrainMtrx, give1dStressMtrx are
 *   introduced since form of this matrices is well known, and for
 *   faster response mainly in linear elastic problems.
 * - Returning a material property (method 'give'). Only for non-standard elements.
 * - Returning real stress state vector(tensor) at gauss point for 3d - case.
 */
class MisesMatTension : public StructuralMaterial, public QCMaterialExtensionInterface
{
protected:
    /// Reference to the basic elastic material.
    IsotropicLinearElasticMaterial linearElasticMaterial;

    /// Elastic shear modulus.
    double G = 0.;

    /// Elastic bulk modulus.
    double K = 0.;

    /// Hardening modulus.
    double H = 0.;

    /// Initial (uniaxial) yield stress.
    ScalarFunction sig0;
    double s0 = 0.;///

    /// critical(maximal) damage.
    double omega_crit = 0.;
    /// exponent in damage function.
    double a = 0.;

    /// tolerance for the yield function in RRM algorithm.
    double yieldTol = 0.;


public:
    MisesMatTension(int n, Domain *d);

    void performPlasticityReturn(const FloatArray &totalStrain, GaussPoint *gp, TimeStep *tStep) const;
    void performPlasticityReturn_PlaneStress(const FloatArrayF<3> &totalStrain, GaussPoint *gp, TimeStep *tStep) const;

    double computeYieldStress(double kappa, GaussPoint *gp, TimeStep *tStep) const;
    double computeYieldStressPrime(double kappa) const;

    double computeDamage(GaussPoint *gp, TimeStep *tStep) const;
    double computeDamageParam(double tempKappa) const;
    double computeDamageParamPrime(double tempKappa) const;
    virtual double computeCumPlastStrain(GaussPoint *gp, TimeStep *tStep) const;

    void initializeFrom(InputRecord &ir) override;

    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) const override { return false; }

    const char *giveInputRecordName() const override { return _IFT_MisesMatTension_Name; }
    const char *giveClassName() const override { return "MisesMatTension"; }

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

    FloatMatrixF<6,6> give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    FloatMatrixF<3,3> givePlaneStressStiffMtrx(MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep) const override;

    FloatMatrixF<1,1> give1dStressStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    FloatArrayF<6> giveRealStressVector_3d(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const override;

    FloatArrayF<3> giveRealStressVector_PlaneStress(const FloatArrayF<3> &totalStrain, GaussPoint *gp,TimeStep *tStep) const override;

    FloatArrayF<1> giveRealStressVector_1d(const FloatArrayF<1> &reducedE, GaussPoint *gp, TimeStep *tStep) const override;

    double give(int aProperty, GaussPoint *gp, TimeStep *tStep) const;
    double giveTemperature(GaussPoint *gp, TimeStep *tStep) const;

    Interface *giveInterface(InterfaceType t) override {
        if ( t == QCMaterialExtensionInterfaceType ) {
            return static_cast< QCMaterialExtensionInterface * >(this);
        } else {
            return nullptr;
        }
    }
    double giveQcElasticParamneter() override { return 9*K*G/(3*K+G); } // E
    double giveQcPlasticParamneter() override { return s0; } // sig0
    double giveQcPlasticHardeningParamneter() override { return H; }
    void computeStatusDataFromStrain(FloatArray &statusData, FloatArray strain) override;
    double giveValueOfQcAdapriveRefinementCriterion(Element *e, int criterionType) override { return -1.; } // adaptive criterion is not implemented
    

protected:
    void computeGLPlasticStrain(const FloatMatrix &F, FloatMatrix &Ep, FloatMatrix b, double J);


    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;
};

//=============================================================================


class MisesMatTensionStatus : public StructuralMaterialStatus
{
protected:
    /// Plastic strain (initial).
    FloatArray plasticStrain;

    /// Plastic strain (final).
    FloatArray tempPlasticStrain;

    /// Deviatoric trial stress - needed for tangent stiffness.
    FloatArray trialStressD;

    /// volumetric trial stress - needed for tangent stiffness.
    double trialStressV = 0.;

    FloatArray effStress;
    FloatArray tempEffStress;

    /// Cumulative plastic strain (initial).
    double kappa = 0.;

    /// Cumulative plastic strain (final).
    double tempKappa = 0.;

    /// damage variable (initial).
    double damage = 0.;

    /// damage variable (final).
    double tempDamage = 0.;


public:
    MisesMatTensionStatus(GaussPoint *g);

    const FloatArray &givePlasticStrain() const { return plasticStrain; }

    const FloatArray &giveTrialStressDev() const { return trialStressD; }

    double giveTrialStressVol() const { return trialStressV; }

    double giveDamage() const { return damage; }
    double giveTempDamage() const { return tempDamage; }

    double giveCumulativePlasticStrain() const { return kappa; }
    double giveTempCumulativePlasticStrain() const { return tempKappa; }

    const FloatArray &giveTempEffectiveStress() const { return tempEffStress; }
    const FloatArray &giveEffectiveStress() const { return effStress; }

    void letTempPlasticStrainBe(const FloatArray &values) { tempPlasticStrain = values; }
    const FloatArray &getTempPlasticStrain() const { return tempPlasticStrain; }

    void letTrialStressDevBe(const FloatArray &values) { trialStressD = values; }

    void letEffectiveStressBe(const FloatArray &values) { effStress = values; }

    void letTempEffectiveStressBe(const FloatArray &values) { tempEffStress = values; }


    void setTrialStressVol(double value) { trialStressV = value; }


    void setTempCumulativePlasticStrain(double value) { tempKappa = value; }

    void setTempDamage(double value) { tempDamage = value; }

    const FloatArray &givePlasDef() { return plasticStrain; }

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    void initTempStatusFromData(const FloatArray &statusData) override;
    void initTempStatus() override;

    void updateYourself(TimeStep *tStep) override;

    const char *giveClassName() const override { return "MisesMatTensionStatus"; }


    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;
};
} // end namespace oofem
#endif // misesmattension_h
