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
 *               Copyright (C) 1993 - 2014   Borek Patzak
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

#ifndef microplane_m1_aniso_h
#define microplane_m1_aniso_h

// by commenting this out we can switch back to the old implementation (2D version not inheriting from MicroplaneMaterial)
//#define microplane_m1_aniso_new_implementation

#ifdef microplane_m1_aniso_new_implementation
// ========================= new implementation =========================

#include "structuralms.h"
#include "microplanematerial.h"
#include "gausspoint.h"


///@name Input fields for M1anisoMaterial
//@{
 #define _IFT_M1anisoMaterial_Name "microplane_m1_aniso"
 #define _IFT_M1anisoMaterial_s0 "s0"
 #define _IFT_M1anisoMaterial_hn "hn"
 #define _IFT_M1anisoMaterial_talpha "talpha"
//@}

namespace oofem {
class M1anisoMaterialStatus : public StructuralMaterialStatus
{
protected:
    FloatArray sigN, tempSigN, epspN, tempEpspN;
    IntArray plasticState;

public:
    M1anisoMaterialStatus(GaussPoint *g);

    const char *giveClassName() const override { return "M1anisoMaterialStatus"; }
    void letTempNormalMplaneStressesBe(FloatArray sigmaN) { tempSigN =  std :: move(sigmaN); }
    void letTempNormalMplanePlasticStrainsBe(FloatArray epsilonpN) { tempEpspN =  std :: move(epsilonpN); }
    void letPlasticStateIndicatorsBe(IntArray plSt) { plasticState =  std :: move(plSt); }
    const FloatArray &giveNormalMplaneStresses() { return sigN; }
    const FloatArray &giveTempNormalMplaneStresses() { return tempSigN; }
    const FloatArray &giveNormalMplanePlasticStrains() { return epspN; }
    const FloatArray &giveTempNormalMplanePlasticStrains() { return tempEpspN; }
    const IntArray &givePlasticStateIndicators() { return plasticState; }
    void initTempStatus() override;
    void printOutputAt(FILE *file, TimeStep *tStep) const override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

};

/**
 * Simple microplane model - version M1, just with normal microplane strains.
 */
class M1anisoMaterial : public MicroplaneMaterial
{
protected:
    double EN = 0.; // normal microplane elastic modulus
    double ENtan = 0.; // normal microplane tangent (elastoplastic) modulus
    double HN = 0.; // normal microplane hardening/softening modulus
    double s0 = 0.; // normal microplane initial yield stress

public:
    /**
     * Constructor. Creates M1anisoMaterial belonging to domain d, with number n.
     * @param n Material number.
     * @param d Domain to which newly created material belongs.
     */
    M1anisoMaterial(int n, Domain *d);

    FloatArrayF<6> giveRealStressVector_3d(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<6,6> give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    const char *giveClassName() const override { return "M1anisoMaterial"; }
    const char *giveInputRecordName() const override { return _IFT_M1anisoMaterial_Name; }
    void initializeFrom(InputRecord &ir) override;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override { return new M1anisoMaterialStatus(gp); }

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;
};
} // end namespace oofem

#else
// ========================= old implementation =========================

 #include "structuralms.h"
 #include "structuralmaterial.h"
 #include "gausspoint.h"

 #include "sm/Materials/qcmaterialextensioninterface.h"


///@name Input fields for M1anisoMaterial
//@{
 #define _IFT_M1anisoMaterial_Name "microplane_m1_aniso"
 #define _IFT_M1anisoMaterial_men "men"
 #define _IFT_M1anisoMaterial_ms0 "ms0"
 #define _IFT_M1anisoMaterial_mhn "mhn"
 #define _IFT_M1anisoMaterial_talpha "talpha"
 #define _IFT_M1anisoMaterial_ma "ma"
 #define _IFT_M1anisoMaterial_mw "mw"
//@}

namespace oofem {
/**
 * Status of simple microplane model - version M1, just with normal microplane strains.
 */
class M1anisoMaterialStatus : public StructuralMaterialStatus
{
protected:
    FloatArray sigN, tempSigN, sigNyield;

public:
    M1anisoMaterialStatus(GaussPoint *g);
    virtual ~M1anisoMaterialStatus();

    // definition
    const char *giveClassName() const override { return "M1anisoMaterialStatus"; }
    void letTempNormalMplaneStressesBe(FloatArray sigmaN) { tempSigN =  std :: move(sigmaN); }
    void letNormalMplaneYieldStressesBe(FloatArray sigmaNyield) { sigNyield =  std :: move(sigmaNyield); }
    const FloatArray &giveNormalMplaneStresses() { return sigN; }
    const FloatArray &giveTempNormalMplaneStresses() { return tempSigN; }
    const FloatArray &giveNormalMplaneYieldStresses() { return sigNyield; }
    void initTempStatus() override;
    void printOutputAt(FILE *file, TimeStep *tStep) const override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;
};


/**
 * Simple microplane model - version M1, just with normal microplane strains.
 */
 class M1anisoMaterial : public StructuralMaterial, public QCMaterialExtensionInterface
{
protected:
    //double E = 0.; // Young's modulus
    double nu = 0.; // Poisson ratio
    //double EN = 0.; // normal microplane elastic modulus
    //double HN = 0.; // normal microplane hardening/softening modulus
    //double s0 = 0.; // normal microplane initial yield stress
    int nmp; // number of microplanes
    FloatMatrix n; // microplane normals
    FloatMatrix N; // N = n x n in Voigt notation
    FloatMatrix NN; // NN = n x n x n x n in special notation
    FloatArray mw; // microplane weights
    FloatArray ma; // microplane angles
    FloatArray mEN; // microplane normal microplane elastic modulus
    FloatArray HN; // normal microplane hardening/softening modulus
    FloatArray s0; // normal microplane initial yield stress
public:
    /**
     * Constructor. Microplane Material M1 belonging to domain d, with number n.
     * @param n Material number.
     * @param d Domain to which newly created material belongs.
     */
    M1anisoMaterial(int n, Domain *d);

    //void giveRealStressVector_PlaneStress(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep);
    FloatArrayF< 3 > giveRealStressVector_PlaneStress(const FloatArrayF< 3 > &totalStrain, GaussPoint *gp, TimeStep *tStep) const override;

    FloatMatrixF<3,3> giveElasticPlaneStressStiffMtrx() const;
    //void givePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);
    FloatMatrixF<3,3> givePlaneStressStiffMtrx(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;

    const char *giveClassName() const override { return "M1anisoMaterial"; }
    const char *giveInputRecordName() const override { return _IFT_M1anisoMaterial_Name; }
    void initializeFrom(InputRecord &ir) override;
    bool hasMaterialModeCapability(MaterialMode mode) const override;
    MaterialStatus *CreateStatus(GaussPoint *gp) const override { return new M1anisoMaterialStatus(gp); }

    double giveQcElasticParamneter() override { return std :: numeric_limits< float > :: infinity(); }
    double giveQcPlasticParamneter() override { return std :: numeric_limits< float > :: infinity(); }
    double giveQcPlasticHardeningParamneter() override { return std :: numeric_limits< float > :: infinity(); }
    void computeStatusDataFromStrain(FloatArray &statusData, FloatArray strain) override {OOFEM_ERROR("This material is not implemented to be homogenized")};
    double giveValueOfQcAdapriveRefinementCriterion(Element *e, int criterionType) override;

    Interface *giveInterface(InterfaceType t) override {
      if ( t == QCMaterialExtensionInterfaceType ) {
	return static_cast< QCMaterialExtensionInterface * >(this);
      } else {
	return nullptr;
      }
    }


protected:
    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

};
} // end namespace oofem
#endif // end of old implementation
#endif // microplane_m1_aniso_h
