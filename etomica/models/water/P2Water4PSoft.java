
package etomica.models.water;

import etomica.api.IAtomPositioned;
import etomica.api.IAtomSet;
import etomica.api.IVector;
import etomica.api.IVector3D;
import etomica.atom.MoleculeOrientedDynamic;
import etomica.potential.IPotentialTorque;
import etomica.space.ISpace;
import etomica.space.Tensor;

/** 
 * 3-point potential for water that can calculate gradient and torque (for the
 * center of mass of the water molecule).
 */
public class P2Water4PSoft extends P2Water4P implements IPotentialTorque {

    public P2Water4PSoft(ISpace space, double sigma, double epsilon,
            double chargeM, double chargeH) {
        super(space, sigma, epsilon, chargeM, chargeH);
        gradient = new IVector3D[2];
        gradient[0] = (IVector3D)space.makeVector();
        gradient[1] = (IVector3D)space.makeVector();
		torque = new IVector3D[2];
		torque[0] = (IVector3D)space.makeVector();
		torque[1] = (IVector3D)space.makeVector();
		fWork = (IVector3D)space.makeVector();
		gradientAndTorque = new IVector[][]{gradient,torque};
		epsilon48 = epsilon*48.0;
    }

    public IVector[][] gradientAndTorque(IAtomSet pair){
        MoleculeOrientedDynamic water1 = (MoleculeOrientedDynamic)pair.getAtom(0);
        MoleculeOrientedDynamic water2 = (MoleculeOrientedDynamic)pair.getAtom(1);

        //compute O-O distance to consider truncation	
        IVector O1r = ((IAtomPositioned)water1.getChildList().getAtom(2)).getPosition();
        IVector O2r = ((IAtomPositioned)water2.getChildList().getAtom(2)).getPosition();

        work.Ev1Mv2(O1r, O2r);
        shift.Ea1Tv1(-1,work);
        boundary.nearestImage(work);
        shift.PE(work);
        double r2 = work.squared();

        if(r2<1.6) {
            throw new RuntimeException("waters are overlapped in gradient");
        }

        double s2 = sigma2/r2;
        double s6 = s2*s2*s2;
        double du = -epsilon48*s6*(s6 - 0.5);
	
        gradient[0].Ea1Tv1(du/r2,work);

        IVector com1 = water1.getPosition();
        IVector com2 = water2.getPosition();

        work.Ev1Mv2(O2r, com2);
        work.XE(gradient[0]);
        torque[1].E(work);
        work.Ev1Mv2(com1, O1r);
        work.XE(gradient[0]);
        torque[0].E(work);

		IVector H11r = ((IAtomPositioned)water1.getChildList().getAtom(0)).getPosition();
		IVector H12r = ((IAtomPositioned)water1.getChildList().getAtom(1)).getPosition();
		IVector H21r = ((IAtomPositioned)water2.getChildList().getAtom(0)).getPosition();
		IVector H22r = ((IAtomPositioned)water2.getChildList().getAtom(1)).getPosition();
        IVector M1r = ((IAtomPositioned)water1.getChildList().getAtom(3)).getPosition();
        IVector M2r = ((IAtomPositioned)water2.getChildList().getAtom(3)).getPosition();

        // M1-M2
        work.Ev1Mv2(M1r, M2r);
        work.PE(shift);
        r2 = work.squared();
        fWork.Ea1Tv1(-chargeMM/(r2*Math.sqrt(r2)), work);
        gradient[0].PE(fWork);
        work.Ev1Mv2(M2r, com2);
        work.XE(fWork);
        torque[1].PE(work);
        work.Ev1Mv2(com1, M1r);
        work.XE(fWork);
        torque[0].PE(work);

        // M1-H21
        work.Ev1Mv2(M1r, H21r);
        work.PE(shift);
        r2 = work.squared();
        fWork.Ea1Tv1(-chargeMH/(r2*Math.sqrt(r2)), work);
        gradient[0].PE(fWork);
        work.Ev1Mv2(H21r, com2);
        work.XE(fWork);
        torque[1].PE(work);
        work.Ev1Mv2(com1, M1r);
        work.XE(fWork);
        torque[0].PE(work);

        // M1-H22
        work.Ev1Mv2(M1r, H22r);
        work.PE(shift);
        r2 = work.squared();
        fWork.Ea1Tv1(-chargeMH/(r2*Math.sqrt(r2)), work);
        gradient[0].PE(fWork);
        work.Ev1Mv2(H22r, com2);
        work.XE(fWork);
        torque[1].PE(work);
        work.Ev1Mv2(com1, M1r);
        work.XE(fWork);
        torque[0].PE(work);

        // H11-M2
        work.Ev1Mv2(H11r, M2r);
        work.PE(shift);
        r2 = work.squared();
        fWork.Ea1Tv1(-chargeMH/(r2*Math.sqrt(r2)), work);
        gradient[0].PE(fWork);
        work.Ev1Mv2(M2r, com2);
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

        // H12-M2
        work.Ev1Mv2(H12r, M2r);
        work.PE(shift);
        r2 = work.squared();
        fWork.Ea1Tv1(-chargeMH/(r2*Math.sqrt(r2)), work);
        gradient[0].PE(fWork);
        work.Ev1Mv2(M2r, com2);
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
    
    public IVector[] gradient(IAtomSet atoms) {
        // do extra work to calculate torque
        gradientAndTorque(atoms);
        return gradient;
    }
    
    public IVector[] gradient(IAtomSet atoms, Tensor pressureTensor) {
        gradientAndTorque(atoms);
        //FIXME
        //pressureTensor.PEv1v2(gradient[0],dr);
        return gradient;
    }

    public double virial(IAtomSet atoms) {
        //FIXME
        return 0;
    }

    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }
    
	public double getSigma() {return sigma;}
    
	public double getEpsilon() {return epsilon;}
	
    private static final long serialVersionUID = 1L;
	protected final IVector3D[] gradient, torque;
	protected final IVector[][] gradientAndTorque;
	protected double epsilon48;
	protected final IVector3D fWork;
}
