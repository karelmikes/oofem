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
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef hemokunzelmat_h
#define hemokunzelmat_h

#include "tm/Materials/transportmaterial.h"

///@name Input fields for HeMoKunzelMaterial
//@{
#define _IFT_HeMoKunzelMaterial_Name "hemokunzel"
#define _IFT_HeMoKunzelMaterial_iso_type "iso_type"
#define _IFT_HeMoKunzelMaterial_permeability_type "permeability_type"
#define _IFT_HeMoKunzelMaterial_iso_a "iso_a"
#define _IFT_HeMoKunzelMaterial_iso_n "iso_n"
#define _IFT_HeMoKunzelMaterial_iso_wh "iso_wh"
#define _IFT_HeMoKunzelMaterial_iso_b "iso_b"
#define _IFT_HeMoKunzelMaterial_mu "mu"
#define _IFT_HeMoKunzelMaterial_pl "pl"
#define _IFT_HeMoKunzelMaterial_a "a"
#define _IFT_HeMoKunzelMaterial_lambda0 "lambda0"
#define _IFT_HeMoKunzelMaterial_b "b"
#define _IFT_HeMoKunzelMaterial_rhoh2o "rhoh2o"
#define _IFT_HeMoKunzelMaterial_cs "cs"
#define _IFT_HeMoKunzelMaterial_cw "cw"
#define _IFT_HeMoKunzelMaterial_hv "hv"
#define _IFT_HeMoKunzelMaterial_perm_h "perm_h"
#define _IFT_HeMoKunzelMaterial_perm_dwh "perm_dw(h)"
#define _IFT_HeMoKunzelMaterial_perm_wv "perm_wv"
#define _IFT_HeMoKunzelMaterial_perm_dwwv "perm_dw(wv)"
//@}

namespace oofem {
/**
 */
class HeMoKunzelMaterial : public TransportMaterial
{
protected:

    enum isothermType { Hansen, Kunzeliso } Isotherm;

    enum permeabilityType { Multilin_h, Multilin_wV, Kunzelperm } Permeability;
    /// values of the multilinear permeability
    FloatArray perm_h;
    FloatArray perm_Dwh;
    FloatArray perm_wV;
    FloatArray perm_DwwV;

    double A;            ///< water absorption coefficient [kg m^-2 s^-0.5]

    double iso_wh;       ///< Parameter of Hansen's/Kunzel's isotherm - max. adsorbed water content [kg/m^3]
    double iso_n;        ///< Parameter of Hansen's isotherm - exponent
    double iso_a;        ///< Parameter of Hansen's isotherm
    double iso_b;        ///< Parameter of Kunzel's isotherm

    double mu;           ///< water vapor diffusion resistance [-]
    double PL;           ///< ambient atmospheric pressure [Pa]

    double lambda0;      ///< thermal conductivity [W m^-1 K^-1]
    double b;            ///< thermal conductivity supplement [-]
    double rho;          ///< bulk density of dry building material [kg m^-3]
    double rhoH2O;       ///< water density [kg m^-3]
    double cs;           ///< specific heat capacity of the building material [J kg^-1 K^-1]
    double cw;           ///< specific heat capacity of liquid water [J kg^-1 K^-1]
    double hv;           ///< latent heat of phase change/evaporation enthalpy of pure water [J/kg]

public:
    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n Material number.
     * @param d Domain to which new material will belong.
     */
    HeMoKunzelMaterial(int n, Domain *d) : TransportMaterial(n, d) { }

    void giveFluxVector(FloatArray &answer, GaussPoint *gp, const FloatArray &grad, const FloatArray &field, TimeStep *tStep) override;

    void giveCharacteristicMatrix(FloatMatrix &answer,
                                  MatResponseMode mode,
                                  GaussPoint *gp,
                                  TimeStep *atTime) override;

    double giveCharacteristicValue(MatResponseMode mode,
                                   GaussPoint *gp,
                                   TimeStep *atTime) override;

    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) const override;

    IRResultType initializeFrom(InputRecord *ir) override;

    double give(int aProperty, GaussPoint *gp) const override;

    bool hasMaterialModeCapability(MaterialMode mode) const override;

    const char *giveInputRecordName() const override { return _IFT_HeMoKunzelMaterial_Name; }
    const char *giveClassName() const override { return "HeMoKunzelMaterial"; }

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *atTime) override;
    double giveHumidity(GaussPoint *gp, ValueModeType mode) const override;

    /// returns water content (in kg/m^3)
    double giveMoistureContent(double h) const;

    /// computes derivative of the moisture storage function (sorption isotherm) with respect to relative humidity
    double giveMoistureContentDerivative(double h) const;

    double computeWaterVaporPerm(double T) const;
    double computeSatVaporPressure(double T) const;
    double computeSatVaporPressureDerivative(double T) const;
    double computeDw(double h) const;

protected:
    void computeConductivityMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime) const;
    void matcond1d(FloatMatrix &d, GaussPoint *gp, MatResponseMode mode, TimeStep *atTime) const;
    void matcond2d(FloatMatrix &d, GaussPoint *gp, MatResponseMode mode, TimeStep *atTime) const;
    void matcond3d(FloatMatrix &d, GaussPoint *gp, MatResponseMode mode, TimeStep *atTime) const;

    double computeCapacityCoeff(MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);

    double perm_mm(double h, double T) const;
    double perm_mh(double h, double T) const;
    double perm_hm(double h, double T) const;
    double perm_hh(double h, double T) const;

};
} // end namespace oofem
#endif // hemokunzelmat_h
