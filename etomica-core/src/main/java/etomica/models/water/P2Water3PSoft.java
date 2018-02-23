/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


package etomica.models.water;

import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeOrientedDynamic;
import etomica.potential.IPotentialMolecularTorque;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;

/** 
 * 3-point potential for water that can calculate gradient and torque (for the
 * center of mass of the water molecule).
 */
public class P2Water3PSoft extends P2Water3P implements IPotentialMolecularTorque {

	public P2Water3PSoft(Space space, double sigma, double epsilon,
                         double chargeO, double chargeH) {
		super(space, sigma, epsilon, chargeO, chargeH);
		gradient = new Vector[2];
		gradient[0] = space.makeVector();
		gradient[1] = space.makeVector();
        torque = new Vector[2];
        torque[0] = space.makeVector();
        torque[1] = space.makeVector();
        fWork = space.makeVector();
        gradientAndTorque = new Vector[][]{gradient,torque};
        epsilon48 = epsilon*48.0;
	}

    public Vector[][] gradientAndTorque(IMoleculeList pair){
		MoleculeOrientedDynamic water1 = (MoleculeOrientedDynamic)pair.get(0);
		MoleculeOrientedDynamic water2 = (MoleculeOrientedDynamic)pair.get(1);
		
		//compute O-O distance to consider truncation	
		Vector O1r = (water1.getChildList().get(2)).getPosition();
		Vector O2r = (water2.getChildList().get(2)).getPosition();

		work.Ev1Mv2(O1r, O2r);
        shift.Ea1Tv1(-1,work);
		boundary.nearestImage(work);
        shift.PE(work);
		double r2 = work.squared();

		if(r2<1.6) {
		    throw new RuntimeException("waters are overlapped in gradient");
		}
	
		double du = -chargeOO/Math.sqrt(r2);

        double s2 = sigma2/r2;
        double s6 = s2*s2*s2;
        du += -epsilon48*s6*(s6 - 0.5);
	
        gradient[0].Ea1Tv1(du/r2,work);

        Vector com1 = water1.getPosition();
        Vector com2 = water2.getPosition();

        work.Ev1Mv2(O2r, com2);
        work.XE(gradient[0]);
        torque[1].E(work);
        
		Vector H11r = water1.getChildList().get(0).getPosition();
		Vector H12r = water1.getChildList().get(1).getPosition();
		Vector H21r = water2.getChildList().get(0).getPosition();
		Vector H22r = water2.getChildList().get(1).getPosition();

        // O1-H21
        work.Ev1Mv2(O1r, H21r);
        work.PE(shift);
        r2 = work.squared();
        fWork.Ea1Tv1(-chargeOH/(r2*Math.sqrt(r2)), work);
        gradient[0].PE(fWork);
        work.Ev1Mv2(H21r, com2);
        work.XE(fWork);
        torque[1].PE(work);

        // O1-H22
        work.Ev1Mv2(O1r, H22r);
        work.PE(shift);
        r2 = work.squared();
        fWork.Ea1Tv1(-chargeOH/(r2*Math.sqrt(r2)), work);
        work.Ev1Mv2(H22r, com2);
        work.XE(fWork);
        torque[1].PE(work);

        gradient[0].PE(fWork);

        // calculate torque for O1 (f on O1 is -gradient[0])
        work.Ev1Mv2(com1, O1r);
        work.XE(gradient[0]);
        torque[0].E(work);

        // H11-O2
        work.Ev1Mv2(H11r, O2r);
        work.PE(shift);
        r2 = work.squared();
        fWork.Ea1Tv1(-chargeOH/(r2*Math.sqrt(r2)), work);
        gradient[0].PE(fWork);
        work.Ev1Mv2(O2r, com2);
        work.XE(fWork);
        torque[1].PE(work);
        work.Ev1Mv2(com1, H11r);
        work.XE(fWork);
        torque[0].PE(work);

        // H11-H21
        work.Ev1Mv2(H11r, H21r);
        work.PE(shift);
        r2 = work.squared();
        fWork.Ea1Tv1(-chargeHH/(r2*Math.sqrt(r2)), work);
        gradient[0].PE(fWork);
        work.Ev1Mv2(H21r, com2);
        work.XE(fWork);
        torque[1].PE(work);
        work.Ev1Mv2(com1, H11r);
        work.XE(fWork);
        torque[0].PE(work);

        // H11-H22
        work.Ev1Mv2(H11r, H22r);
        work.PE(shift);
        r2 = work.squared();
        fWork.Ea1Tv1(-chargeHH/(r2*Math.sqrt(r2)), work);
        gradient[0].PE(fWork);
        work.Ev1Mv2(H22r, com2);
        work.XE(fWork);
        torque[1].PE(work);
        work.Ev1Mv2(com1, H11r);
        work.XE(fWork);
        torque[0].PE(work);

        // H12-O2
        work.Ev1Mv2(H12r, O2r);
        work.PE(shift);
        r2 = work.squared();
        fWork.Ea1Tv1(-chargeOH/(r2*Math.sqrt(r2)), work);
        gradient[0].PE(fWork);
        work.Ev1Mv2(O2r, com2);
        work.XE(fWork);
        torque[1].PE(work);
        work.Ev1Mv2(com1, H12r);
        work.XE(fWork);
        torque[0].PE(work);

        // H12-H21
        work.Ev1Mv2(H12r, H21r);
        work.PE(shift);
        r2 = work.squared();
        fWork.Ea1Tv1(-chargeHH/(r2*Math.sqrt(r2)), work);
        gradient[0].PE(fWork);
        work.Ev1Mv2(H21r, com2);
        work.XE(fWork);
        torque[1].PE(work);
        work.Ev1Mv2(com1, H12r);
        work.XE(fWork);
        torque[0].PE(work);

        // H12-H22
        work.Ev1Mv2(H12r, H22r);
        work.PE(shift);
        r2 = work.squared();
        fWork.Ea1Tv1(-chargeHH/(r2*Math.sqrt(r2)), work);
        gradient[0].PE(fWork);
        work.Ev1Mv2(H22r, com2);
        work.XE(fWork);
        torque[1].PE(work);
        work.Ev1Mv2(com1, H12r);
        work.XE(fWork);
        torque[0].PE(work);

        gradient[1].Ea1Tv1(-1, gradient[0]);

		return gradientAndTorque;
	}
    
    public Vector[] gradient(IMoleculeList atoms) {
        // do extra work to calculate torque
        gradientAndTorque(atoms);
        return gradient;
    }
    
    public Vector[] gradient(IMoleculeList atoms, Tensor pressureTensor) {
        gradientAndTorque(atoms);
        //FIXME
        //pressureTensor.PEv1v2(gradient[0],dr);
        return gradient;
    }

    public double virial(IMoleculeList atoms) {
        //FIXME
        return 0;
    }

    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }
    
	public double getSigma() {return sigma;}
    
	public double getEpsilon() {return epsilon;}
	
    private static final long serialVersionUID = 1L;
	protected final Vector[] gradient, torque;
	protected final Vector[][] gradientAndTorque;
	protected double epsilon48;
	protected final Vector fWork;
}
