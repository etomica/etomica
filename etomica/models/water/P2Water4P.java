/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


package etomica.models.water;

import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.IVectorMutable;
import etomica.potential.PotentialMolecular;
import etomica.space.ISpace;

/** 
 * 4-point potential for water.  Potential parameters are typically defined
 * by a convenience subclass.
 * 
 * @author David Kofke, Andrew Schultz
 */
public class P2Water4P extends PotentialMolecular {

	public P2Water4P(ISpace space, double sigma, double epsilon, double chargeM, double chargeH) {
		super(2, space);
        this.sigma = sigma;
        sigma2 = sigma*sigma;
		work = space.makeVector();
		shift = space.makeVector();
        this.epsilon = epsilon;
        epsilon4 = 4*epsilon;
        chargeMM = chargeM * chargeM;
        chargeMH = chargeM * chargeH;
        chargeHH = chargeH * chargeH;
        this.chargeM = chargeM;
        this.chargeH = chargeH;
	}

    public void setBox(IBox box) {
        boundary = box.getBoundary();
    }

    public double energy(IMoleculeList pair){
		double sum = 0.0;
		double r2 = 0.0;

		IMolecule water1 = pair.getMolecule(0);
		IMolecule water2 = pair.getMolecule(1);
		
		//compute O-O distance to consider truncation	
        IVectorMutable O1r = (water1.getChildList().getAtom(2)).getPosition();
        IVectorMutable O2r = (water2.getChildList().getAtom(2)).getPosition();

		work.Ev1Mv2(O1r, O2r);
        shift.Ea1Tv1(-1,work);
		boundary.nearestImage(work);
        shift.PE(work);
        final boolean zeroShift = shift.squared() < 0.1;
		r2 = work.squared();

		if(r2<1.6) return Double.POSITIVE_INFINITY;
	
		double s2 = sigma2/(r2);
		double s6 = s2*s2*s2;
		sum += epsilon4*s6*(s6 - 1.0);
		
        IVectorMutable H11r = water1.getChildList().getAtom(0).getPosition();
        IVectorMutable H12r = water1.getChildList().getAtom(1).getPosition();
        IVectorMutable H21r = water2.getChildList().getAtom(0).getPosition();
        IVectorMutable H22r = water2.getChildList().getAtom(1).getPosition();
        IVectorMutable M1r = water1.getChildList().getAtom(3).getPosition();
        IVectorMutable M2r = water2.getChildList().getAtom(3).getPosition();
        		
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
	protected IBoundary boundary;
	protected final double chargeH;
	protected final double chargeM;
	protected final double chargeMM, chargeMH, chargeHH;
	protected final IVectorMutable work, shift;
}
