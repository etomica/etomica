/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


package etomica.models.water;

import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.IMoleculePositionDefinition;
import etomica.potential.PotentialMolecular;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;

/** 
 * 4-point potential for water.  Potential parameters are typically defined
 * by a convenience subclass.
 * 
 * @author David Kofke, Andrew Schultz
 */
public class P2Water4P extends PotentialMolecular {

	public P2Water4P(Space space, double sigma, double epsilon, double chargeH, double rCut, IMoleculePositionDefinition positionDefinition) {
		super(2, space);
        this.sigma = sigma;
        sigma2 = sigma*sigma;
		work = space.makeVector();
		shift = space.makeVector();
		com1 = space.makeVector();
		com2 = space.makeVector();
        this.epsilon = epsilon;
        epsilon4 = 4*epsilon;
        this.chargeM = -2*chargeH;
        chargeMM = chargeM * chargeM;
        chargeMH = chargeM * chargeH;
        chargeHH = chargeH * chargeH;
        this.chargeH = chargeH;
        this.rCut = rCut;
        this.positionDefinition = positionDefinition;
	}

    public void setBox(Box box) {
        boundary = box.getBoundary();
    }

    public double energy(IMoleculeList pair){

		IMolecule water1 = pair.get(0);
		IMolecule water2 = pair.get(1);
		
		//compute O-O distance to consider truncation	
        Vector O1r = (water1.getChildList().get(2)).getPosition();
        Vector O2r = (water2.getChildList().get(2)).getPosition();

		work.Ev1Mv2(O1r, O2r);
        shift.Ea1Tv1(-1,work);
		boundary.nearestImage(work);
        shift.PE(work);
        final boolean zeroShift = shift.squared() < 0.1;
		double r2 = work.squared();

		if(r2<1.6) return Double.POSITIVE_INFINITY;
		
		com1.E(positionDefinition.position(water1));
		com2.E(positionDefinition.position(water2));
		work.Ev1Mv2(com1, com2);
		work.PE(shift);
		if(work.squared() > rCut*rCut) return 0;
		
		double s2 = sigma2/(r2);
		double s6 = s2*s2*s2;
		double sum = epsilon4*s6*(s6 - 1.0);
		
        Vector H11r = water1.getChildList().get(0).getPosition();
        Vector H12r = water1.getChildList().get(1).getPosition();
        Vector H21r = water2.getChildList().get(0).getPosition();
        Vector H22r = water2.getChildList().get(1).getPosition();
        Vector M1r = water1.getChildList().get(3).getPosition();
        Vector M2r = water2.getChildList().get(3).getPosition();
        		
        if (zeroShift) {
            r2 = M1r.Mv1Squared(M2r);
            sum += chargeMM/Math.sqrt(r2);
            r2 = M1r.Mv1Squared(H21r);
            sum += chargeMH/Math.sqrt(r2);
            r2 = M1r.Mv1Squared(H22r);
            sum += chargeMH/Math.sqrt(r2);
            r2 = H11r.Mv1Squared(M2r);
            sum += chargeMH/Math.sqrt(r2);
            r2 = H11r.Mv1Squared(H21r);
            sum += chargeHH/Math.sqrt(r2);
            r2 = H11r.Mv1Squared(H22r);
            sum += chargeHH/Math.sqrt(r2);
            r2 = H12r.Mv1Squared(M2r);
            sum += chargeMH/Math.sqrt(r2);
            r2 = H12r.Mv1Squared(H21r);
            sum += chargeHH/Math.sqrt(r2);
            r2 = H12r.Mv1Squared(H22r);
            sum += chargeHH/Math.sqrt(r2);
        }
        else {
            shift.PE(M1r);
            r2 = M2r.Mv1Squared(shift);
            shift.ME(M1r);
            sum += chargeMM/Math.sqrt(r2);

            shift.PE(M1r);
            r2 = H21r.Mv1Squared(shift);
            shift.ME(M1r);
            sum += chargeMH/Math.sqrt(r2);

            shift.PE(M1r);
            r2 = H22r.Mv1Squared(shift);
            shift.ME(M1r);
            sum += chargeMH/Math.sqrt(r2);

            shift.PE(H11r);
            r2 = M2r.Mv1Squared(shift);
            shift.ME(H11r);
            sum += chargeMH/Math.sqrt(r2);

            shift.PE(H11r);
            r2 = H21r.Mv1Squared(shift);
            shift.ME(H11r);
            sum += chargeHH/Math.sqrt(r2);

            shift.PE(H11r);
            r2 = H22r.Mv1Squared(shift);
            shift.ME(H11r);
            sum += chargeHH/Math.sqrt(r2);

            shift.PE(H12r);
            r2 = M2r.Mv1Squared(shift);
            shift.ME(H12r);
            sum += chargeMH/Math.sqrt(r2);

            shift.PE(H12r);
            r2 = H21r.Mv1Squared(shift);
            shift.ME(H12r);
            sum += chargeHH/Math.sqrt(r2);

            shift.PE(H12r);
            r2 = H22r.Mv1Squared(shift);
            shift.ME(H12r);
            sum += chargeHH/Math.sqrt(r2);
        }

		return sum;																					        
	}

    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

	public double getSigma() {return sigma;}

	public double getEpsilon() {return epsilon;}

    private static final long serialVersionUID = 1L;
	public double sigma , sigma2;
	protected double epsilon, epsilon4;
	protected Boundary boundary;
	protected final double chargeH;
	protected final double chargeM;
	protected final double chargeMM, chargeMH, chargeHH;
	protected final double rCut;
	protected final Vector work, shift,com1,com2;
	protected final IMoleculePositionDefinition positionDefinition;
}
