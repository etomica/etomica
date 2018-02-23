/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


package etomica.potential;

import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;

/** 
 * 3-point potential for CO2.  Includes LJ interactions for each pair of atom types.
 * Molecule is assumed to be rigid.
 * Potential parameters are typically defined by a convenience subclass.
 * 
 * @author Andrew Schultz
 */
public class P2CO2EMP extends PotentialMolecular {

	public P2CO2EMP(Space space, double sigmaC, double sigmaCO, double sigmaO, double epsilonC, double epsilonCO, double epsilonO, double chargeC) {
		super(2, null);
        this.sigmaC = sigmaC;
        sigmaC2 = sigmaC*sigmaC;
        this.sigmaO = sigmaO;
        sigmaO2 = sigmaO*sigmaO;
        this.sigmaCO = sigmaCO;
        sigmaCO2 = sigmaCO*sigmaCO;
		work = space.makeVector();
		shift = space.makeVector();
        this.epsilonC = epsilonC;
        this.epsilonO = epsilonO;
        this.epsilonCO = epsilonCO;
        chargeCC = chargeC * chargeC;
        double chargeO = -0.5*chargeC;
        chargeCO = chargeC * chargeO;
        chargeOO = chargeO * chargeO;
	}

    public void setBox(Box box) {
        boundary = box.getBoundary();
    }

    public double energy(IMoleculeList pair){

		IMolecule m1 = pair.getMolecule(0);
		IMolecule m2 = pair.getMolecule(1);
		
		//compute C-C distance to consider truncation	
        Vector C1r = (m1.getChildList().get(0)).getPosition();
        Vector C2r = (m2.getChildList().get(0)).getPosition();

		work.Ev1Mv2(C1r, C2r);
        shift.Ea1Tv1(-1,work);
		boundary.nearestImage(work);
        shift.PE(work);
        final boolean zeroShift = shift.squared() < 0.1;
		double r2 = work.squared();

		if(r2<1.6) return Double.POSITIVE_INFINITY;
	
		double sum = chargeCC/Math.sqrt(r2);
		double s2 = sigmaC2/r2;
		double s6 = s2*s2*s2;
		sum += 4*epsilonC*s6*(s6 - 1.0);
		
        Vector O11r = (m1.getChildList().get(1)).getPosition();
        Vector O12r = (m1.getChildList().get(2)).getPosition();
        Vector O21r = (m2.getChildList().get(1)).getPosition();
        Vector O22r = (m2.getChildList().get(2)).getPosition();

        if (zeroShift) {
            r2 = C1r.Mv1Squared(O21r);
            sum += chargeCO/Math.sqrt(r2);
            s2 = sigmaCO2/r2;
            s6 = s2*s2*s2;
            sum += 4*epsilonCO*s6*(s6 - 1.0);

            r2 = C1r.Mv1Squared(O22r);
            sum += chargeCO/Math.sqrt(r2);
            s2 = sigmaCO2/r2;
            s6 = s2*s2*s2;
            sum += 4*epsilonCO*s6*(s6 - 1.0);

            r2 = O11r.Mv1Squared(C2r);
            sum += chargeCO/Math.sqrt(r2);
            s2 = sigmaCO2/r2;
            s6 = s2*s2*s2;
            sum += 4*epsilonCO*s6*(s6 - 1.0);

            r2 = O11r.Mv1Squared(O21r);
            sum += chargeOO/Math.sqrt(r2);
            s2 = sigmaO2/r2;
            s6 = s2*s2*s2;
            sum += 4*epsilonO*s6*(s6 - 1.0);

            r2 = O11r.Mv1Squared(O22r);
            sum += chargeOO/Math.sqrt(r2);
            s2 = sigmaO2/r2;
            s6 = s2*s2*s2;
            sum += 4*epsilonO*s6*(s6 - 1.0);

            r2 = O12r.Mv1Squared(C2r);
            sum += chargeCO/Math.sqrt(r2);
            s2 = sigmaCO2/r2;
            s6 = s2*s2*s2;
            sum += 4*epsilonCO*s6*(s6 - 1.0);

            r2 = O12r.Mv1Squared(O21r);
            sum += chargeOO/Math.sqrt(r2);
            s2 = sigmaO2/r2;
            s6 = s2*s2*s2;
            sum += 4*epsilonO*s6*(s6 - 1.0);

            r2 = O12r.Mv1Squared(O22r);
            sum += chargeOO/Math.sqrt(r2);
            s2 = sigmaO2/r2;
            s6 = s2*s2*s2;
            sum += 4*epsilonO*s6*(s6 - 1.0);
        }
        else {
            shift.PE(C1r);
            r2 = O21r.Mv1Squared(shift);
            shift.ME(C1r);
            sum += chargeCO/Math.sqrt(r2);
            s2 = sigmaCO2/r2;
            s6 = s2*s2*s2;
            sum += 4*epsilonCO*s6*(s6 - 1.0);

            shift.PE(C1r);
            r2 = O22r.Mv1Squared(shift);
            shift.ME(C1r);
            sum += chargeCO/Math.sqrt(r2);
            s2 = sigmaCO2/r2;
            s6 = s2*s2*s2;
            sum += 4*epsilonCO*s6*(s6 - 1.0);

            shift.PE(O11r);
            r2 = C2r.Mv1Squared(shift);
            shift.ME(O11r);
            sum += chargeCO/Math.sqrt(r2);
            s2 = sigmaCO2/r2;
            s6 = s2*s2*s2;
            sum += 4*epsilonCO*s6*(s6 - 1.0);

            shift.PE(O11r);
            r2 = O21r.Mv1Squared(shift);
            shift.ME(O11r);
            sum += chargeOO/Math.sqrt(r2);
            s2 = sigmaO2/r2;
            s6 = s2*s2*s2;
            sum += 4*epsilonO*s6*(s6 - 1.0);

            shift.PE(O11r);
            r2 = O22r.Mv1Squared(shift);
            shift.ME(O11r);
            sum += chargeOO/Math.sqrt(r2);
            s2 = sigmaO2/r2;
            s6 = s2*s2*s2;
            sum += 4*epsilonO*s6*(s6 - 1.0);

            shift.PE(O12r);
            r2 = C2r.Mv1Squared(shift);
            shift.ME(O12r);
            sum += chargeCO/Math.sqrt(r2);
            s2 = sigmaCO2/r2;
            s6 = s2*s2*s2;
            sum += 4*epsilonCO*s6*(s6 - 1.0);

            shift.PE(O12r);
            r2 = O21r.Mv1Squared(shift);
            shift.ME(O12r);
            sum += chargeOO/Math.sqrt(r2);
            s2 = sigmaO2/r2;
            s6 = s2*s2*s2;
            sum += 4*epsilonO*s6*(s6 - 1.0);

            shift.PE(O12r);
            r2 = O22r.Mv1Squared(shift);
            shift.ME(O12r);
            sum += chargeOO/Math.sqrt(r2);
            s2 = sigmaO2/r2;
            s6 = s2*s2*s2;
            sum += 4*epsilonO*s6*(s6 - 1.0);
        }
		return sum;																					        
	}

    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public double getSigmaC() {
        return sigmaC;
    }

    public double getSigmaO() {
        return sigmaO;
    }

    public double getSigmaCO() {
        return sigmaCO;
    }

    public double getEpsilonC() {
        return epsilonC;
    }

    public double getEpsilonO() {
        return epsilonO;
    }

    public double getEpsilonCO() {
        return epsilonCO;
    }
    
    private static final long serialVersionUID = 1L;
	public double sigmaC, sigmaC2, sigmaO, sigmaO2, sigmaCO, sigmaCO2;
	protected double epsilonC, epsilonO, epsilonCO;
	protected Boundary boundary;
	protected final double chargeCC, chargeCO, chargeOO;
	protected final Vector work, shift;
}
