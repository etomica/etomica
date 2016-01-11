/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.IPotentialMolecular;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.atom.DipoleSource;
import etomica.atom.IAtomOriented;
import etomica.atom.IAtomPositionDefinition;
import etomica.space.ISpace;
import etomica.space.Tensor;

public class P2ReactionFieldDipole extends PotentialMolecular implements PotentialMolecularSoft, IPotentialMolecularSecondDerivative {

    public P2ReactionFieldDipole(ISpace space,IAtomPositionDefinition positionDefinition) {
        super(2, space);
        this.positionDefinition = positionDefinition;
        iDipole = space.makeVector();
        cavityDipole = space.makeVector();
        dr = space.makeVector();
        dr1 = space.makeVector();
        gradientAndTorque = new IVectorMutable[2][2];
        gradientAndTorque[0][0] = space.makeVector();
        gradientAndTorque[0][1] = space.makeVector();
        gradientAndTorque[1][0] = space.makeVector();
        gradientAndTorque[1][1] = space.makeVector();
        secondDerivative = new Tensor [3];
        secondDerivative[0] = space.makeTensor();
		secondDerivative[1] = space.makeTensor();
		secondDerivative[2] = space.makeTensor();
		a = new IVectorMutable[3];
		a[0] = space.makeVector();
		a[1] = space.makeVector();
		a[2] = space.makeVector();
    }

    /**
     * Returns the dipole source used by this object.
     */
    public DipoleSource getDipoleSource() {
        return dipoleSource;
    }

    /**
     * Sets the dipole source used by this object should use.
     */
    public void setDipoleSource(DipoleSource newDipoleSource) {
        dipoleSource = newDipoleSource;
    }
    
    public double getRange() {
        return cutoff;
    }
    
    public void setRange(double newRange) {
        cutoff = newRange;
        cutoff2 = newRange * newRange;
        fac = 2*(epsilon-1)/(2*epsilon+1)/(cutoff2*cutoff);
    }
    
    /**
     * Returns the dielectric constant of the fluid surrounding the cavity.
     */
    public double getDielectric() {
        return epsilon;
    }
    
    /**
     * Sets the dielectric constant of the fluid surrounding the cavity.
     */
    public void setDielectric(double newDielectric) {
        epsilon = newDielectric;
        if (cutoff > 0) {
            fac = 2*(epsilon-1)/(2*epsilon+1)/(cutoff2*cutoff);
//            System.out.println("epsilon = " + epsilon);
//            System.out.println("cutoff = " + cutoff);
//            System.out.println("fac in setDielectric = " +fac);
//            System.exit(2);
//            fac = 1.0;//TODO
        }
    }

    public void setBox(IBox box) {
        boundary = box.getBoundary();
    }

    public double energy(IMoleculeList molecules) {
    	IMolecule molecule0 = molecules.getMolecule(0);
		IMolecule molecule1 = molecules.getMolecule(1);
		dr.E(positionDefinition.position(molecule1));
		dr.ME(positionDefinition.position(molecule0));
		boundary.nearestImage(dr);//rij
		double r2 = dr.squared();
//		System.out.println("r2 = " + r2 + " rCut2 = " + cutoff*cutoff);
		if (r2 > cutoff*cutoff) return 0;
		
        iDipole.E(dipoleSource.getDipole(molecules.getMolecule(0)));
        double idotj = iDipole.dot(dipoleSource.getDipole(molecules.getMolecule(1)));
//        System.out.println("count = " + (++count));
//        System.out.println("iDipole = " + dipoleSource.getDipole(molecules.getMolecule(0)));
//        System.out.println("jDipole = " + dipoleSource.getDipole(molecules.getMolecule(1)));
//        System.out.println("URF = " + -fac*idotj);
//        System.out.println("fac in energy = " + fac);
//        System.out.println("fac = " + 2*(epsilon-1)/(2*epsilon+1)/(cutoff2*cutoff));
//        System.exit(2);
        return -fac*idotj;
    }

    public IVector[][] gradientAndTorque(IMoleculeList molecules) {
//    	System.out.println("fac at begin = " + fac);
    	IMolecule molecule0 = molecules.getMolecule(0);
		IMolecule molecule1 = molecules.getMolecule(1);
		dr.E(positionDefinition.position(molecule1));
		dr.ME(positionDefinition.position(molecule0));
		boundary.nearestImage(dr);
		gradientAndTorque[0][0].E(0);
		gradientAndTorque[0][1].E(0);
		double r2 = dr.squared();
		if (r2 > cutoff*cutoff){
			gradientAndTorque[1][0].E(0);
			gradientAndTorque[1][1].E(0);
			return gradientAndTorque;
		}
		
        iDipole.E(dipoleSource.getDipole(molecule0));
        iDipole.XE(dipoleSource.getDipole(molecule1));
        iDipole.TE(fac);
        gradientAndTorque[1][0].E(iDipole);
        gradientAndTorque[1][1].Ea1Tv1(-1,iDipole);
		
        if(false){//finite difference TODO
        System.out.println("torque =  " + iDipole);
        iDipole.E(dipoleSource.getDipole(molecule0));
        IVectorMutable iDipoleP = space.makeVector();
        IVectorMutable X = space.makeVector();
        IVectorMutable Y = space.makeVector();
        IVectorMutable Z = space.makeVector();
        X.setX(0, 1);
        Y.setX(1, 1);
        Z.setX(2, 1);
        
        double u0 = -fac*iDipole.dot(dipoleSource.getDipole(molecules.getMolecule(1)));
        double prec = 1.0E-4;
        iDipoleP.E(iDipole);
        iDipoleP.PEa1Tv1(-(iDipole.dot(X)), X);
        double norm = Math.sqrt(iDipoleP.squared());
        dr1.E(X);
        dr1.XE(iDipole);
        dr1.normalize();
        dr1.TE(norm*Math.tan(prec));
        iDipoleP.PE(dr1);
        iDipoleP.normalize();
        iDipoleP.TE(norm);
        iDipoleP.PEa1Tv1(iDipole.dot(X), X);
        double ux= -fac*iDipoleP.dot(dipoleSource.getDipole(molecules.getMolecule(1)));
        double x = -(ux-u0)/prec;
        
        
        iDipoleP.E(iDipole);
        iDipoleP.PEa1Tv1(-iDipole.dot(Y), Y);
        norm = Math.sqrt(iDipoleP.squared());
        dr1.E(Y);
        dr1.XE(iDipole);
        dr1.normalize();
        dr1.TE(norm*Math.tan(prec));
        iDipoleP.PE(dr1);
        iDipoleP.normalize();
        iDipoleP.TE(norm);
        iDipoleP.PEa1Tv1(iDipole.dot(Y), Y);
        double uy= -fac*iDipoleP.dot(dipoleSource.getDipole(molecules.getMolecule(1)));
        double y = -(uy-u0)/prec;
        
        iDipoleP.E(iDipole);
        iDipoleP.PEa1Tv1(-iDipole.dot(Z), Z);
        norm = Math.sqrt(iDipoleP.squared());
        dr1.E(Z);
        dr1.XE(iDipole);
        dr1.normalize();
        dr1.TE(norm*Math.tan(prec));
        iDipoleP.PE(dr1);
        iDipoleP.normalize();
        iDipoleP.TE(norm);
        iDipoleP.PEa1Tv1(iDipole.dot(Z), Z);
        double uz= -fac*iDipoleP.dot(dipoleSource.getDipole(molecules.getMolecule(1)));
        double z = -(uz-u0)/prec;
        
        
        System.out.println("torqueP = "+"(" + x + " , " + y + " , " + z + ")");
        System.exit(2);
    }
        return gradientAndTorque;
    }
    
    public Tensor[] secondDerivative(IMoleculeList molecules){
    	IMolecule molecule0 = molecules.getMolecule(0);
		IMolecule molecule1 = molecules.getMolecule(1);
		dr.E(positionDefinition.position(molecule1));
		dr.ME(positionDefinition.position(molecule0));
		boundary.nearestImage(dr);//rij
		double r2 = dr.squared();
		secondDerivative[0].E(0);
		secondDerivative[1].E(0);
		secondDerivative[2].E(0);
		if (r2 > cutoff*cutoff) return secondDerivative;
//		System.out.println("frac in tensor = " + fac);
		iDipole.E(dipoleSource.getDipole(molecule0));
		IVectorMutable jDipole = space.makeVector();
		jDipole.E(dipoleSource.getDipole(molecule1));
		
		double exi = iDipole.getX(0);//ei and ej is the dipole orientation with mu
		double eyi = iDipole.getX(1);
		double ezi = iDipole.getX(2);
		double exj = jDipole.getX(0);
		double eyj = jDipole.getX(1);
		double ezj = jDipole.getX(2);
		
		IVectorMutable deidxi = space.makeVector();
		IVectorMutable deidyi = space.makeVector();
		IVectorMutable deidzi = space.makeVector();
		IVectorMutable dejdxj = space.makeVector();
		IVectorMutable dejdyj = space.makeVector();
		IVectorMutable dejdzj = space.makeVector();
		
		double [] deidxiD = {0,-ezi,eyi};
		double [] deidyiD = {ezi,0,-exi};
		double [] deidziD = {-eyi,exi,0};
		deidxi.E(deidxiD);
		deidyi.E(deidyiD);
		deidzi.E(deidziD);
		
		double [] dejdxjD = {0,-ezj,eyj};
		double [] dejdyjD = {ezj,0,-exj};
		double [] dejdzjD = {-eyj,exj,0};
		dejdxj.E(dejdxjD);
		dejdyj.E(dejdyjD);
		dejdzj.E(dejdzjD);
		
		//debug only
//		System.out.println("r  = " + dr);
//		System.out.println("ei = " + iDipole);
//		System.out.println("ej = " + jDipole);
		
		a[0].E(iDipole);
		a[0].XE(dejdxj);
		a[1].E(iDipole);
		a[1].XE(dejdyj);
		a[2].E(iDipole);
		a[2].XE(dejdzj);
		secondDerivative[0].E(a);
		//TODO
//		System.out.println("ij = \n" + secondDerivative[0]);//debug only
		
		secondDerivative[0].TE(-fac);//ij and ji are the same
		
//		System.out.println("ij = \n" + secondDerivative[0]);//debug only
//		System.exit(2);
//		System.out.println("ijdxj = " + a[0]);
		
		a[0].E(deidxi);
		a[0].XE(jDipole);
		a[1].E(deidyi);
		a[1].XE(jDipole);
		a[2].E(deidzi);
		a[2].XE(jDipole);
		secondDerivative[1].E(a);
		
//		System.out.println("ii = \n" + secondDerivative[1]);//debug only
		
		secondDerivative[1].TE(-fac);//ii
		
		
//		System.out.println("ii = \n" + secondDerivative[1]);//TODO
//		System.out.println("fac = " + fac);
//		System.exit(2);
//		a[0].TE(-fac);
//		System.out.println("iidxi = " + a[0]);
		
		
		a[0].E(dejdxj);
		a[0].XE(iDipole);
		a[1].E(dejdyj);
		a[1].XE(iDipole);
		a[2].E(dejdzj);
		a[2].XE(iDipole);
		secondDerivative[2].E(a);
		
		
//		System.out.println("jj = \n" + secondDerivative[2]); //debug only
		secondDerivative[2].TE(-fac);//jj 
//		System.out.println("-ij*fac = \n" + secondDerivative[0]);//debug only
//		System.out.println("-ii*fac = \n" + secondDerivative[1]);//debug only
//		System.out.println("-jj*fac = \n" + secondDerivative[2]); //debug only
//		System.exit(2);
		
		if(false){//finite difference from dtaudx TODO
			System.out.println("ij = " + secondDerivative[0]);
			System.out.println("ii = " + secondDerivative[1]);
			System.out.println("jj = " + secondDerivative[2]);
			iDipole.E(dipoleSource.getDipole(molecule0));
			jDipole.E(dipoleSource.getDipole(molecule1));
			IVectorMutable iDipoleP = space.makeVector();
			IVectorMutable jDipoleP = space.makeVector();
			IVectorMutable iTorque = space.makeVector();
			IVectorMutable iTorqueP = space.makeVector();
			IVectorMutable X = space.makeVector();
	        IVectorMutable Y = space.makeVector();
	        IVectorMutable Z = space.makeVector();
	        X.setX(0, 1);
	        Y.setX(1, 1);
	        Z.setX(2, 1);
			iDipoleP.E(iDipole);
			jDipoleP.E(jDipole);
			
			iTorque.E(iDipole);
			iTorque.XE(jDipole);
			iTorque.TE(fac);
			
			double prec = 1.0E-5;
			iDipoleP.E(iDipole);
			iDipoleP.PEa1Tv1(-(iDipole.dot(X)), X);
			double norm = Math.sqrt(iDipoleP.squared());
			dr1.E(X);
			dr1.XE(iDipole);
			dr1.normalize();
			dr1.TE(norm*Math.tan(prec));
			iDipoleP.PE(dr1);
			iDipoleP.normalize();
			iDipoleP.TE(norm);
			iDipoleP.PEa1Tv1(iDipole.dot(X), X);
			
			iTorqueP.E(iDipoleP);
			iTorqueP.XE(jDipole);
			iTorqueP.TE(fac);
			
			iTorqueP.ME(iTorque);
			iTorqueP.TE(-1.0/prec);
			System.out.println("diTaudxi = " + iTorqueP);
			
			
			jDipoleP.E(jDipole);
			jDipoleP.PEa1Tv1(-(jDipole.dot(X)), X);
			norm = Math.sqrt(jDipoleP.squared());
			dr1.E(X);
			dr1.XE(jDipole);
			dr1.normalize();
			dr1.TE(norm*Math.tan(prec));
			jDipoleP.PE(dr1);
			jDipoleP.normalize();
			jDipoleP.TE(norm);
			jDipoleP.PEa1Tv1(jDipole.dot(X), X);
			
			iTorqueP.E(iDipole);
			iTorqueP.XE(jDipoleP);
			iTorqueP.TE(fac);
			
			iTorqueP.ME(iTorque);
			iTorqueP.TE(-1.0/prec);
			System.out.println("diTaudxj = " + iTorqueP);
			
			
			IVectorMutable jTorque  = space.makeVector();
			IVectorMutable jTorqueP  = space.makeVector();
			iDipoleP.E(iDipole);
			jDipoleP.E(jDipole);
			
			jTorque.E(jDipole);
			jTorque.XE(iDipole);
			jTorque.TE(fac);
			
			jDipoleP.E(jDipole);
			jDipoleP.PEa1Tv1(-(jDipole.dot(X)), X);
			norm = Math.sqrt(jDipoleP.squared());
			dr1.E(X);
			dr1.XE(jDipole);
			dr1.normalize();
			dr1.TE(norm*Math.tan(prec));
			jDipoleP.PE(dr1);
			jDipoleP.normalize();
			jDipoleP.TE(norm);
			jDipoleP.PEa1Tv1(jDipole.dot(X), X);
			
			jTorqueP.E(jDipoleP);
			jTorqueP.XE(iDipole);
			jTorqueP.TE(fac);
			
			jTorqueP.ME(jTorque);
			jTorqueP.TE(-1.0/prec);
			System.out.println("djTaudxj = " + jTorqueP);
			
			
			
			
			
			System.exit(2);
		}
		
		
    	return secondDerivative; 
    }
    
    

    public IVector[] gradient(IMoleculeList atoms) {
        return gradientAndTorque[0];
    }

    public IVector[] gradient(IMoleculeList atoms, Tensor pressureTensor) {
        return gradient(atoms);
    }

    public double virial(IMoleculeList atoms) {
        return 0;
    }

    /**
     * Returns a 0-body potential that should be added along with this
     * potential.
     */
    public IPotentialMolecular makeP0() {
        return new P0ReactionField(this.space, this);
    }

    protected final IAtomPositionDefinition positionDefinition;
    protected final IVectorMutable iDipole, cavityDipole;
    protected final IVectorMutable dr ,dr1;
    protected DipoleSource dipoleSource;
    protected IBoundary boundary;
    protected double cutoff2, cutoff;
    protected double epsilon;
    protected final IVectorMutable[][] gradientAndTorque;
    protected final Tensor[] secondDerivative;
    protected final IVectorMutable [] a;
    protected double fac;
    protected int count = 0;//TODO
    
    /**
     * A 0-body potential that should be added along with this potential.  The
     * 0-body potential includes the effective self-interaction of the
     * molecules (the molecule induces a dipole in the surrounding fluid, which
     * has an interaction energy with the molecule).  This part of the
     * potential does not result in a gradient or torque on the molecule and is
     * independent of position or orientation.
     */
    public static class P0ReactionField extends PotentialMolecular implements IPotential0Lrc, PotentialMolecularSoft {

        public P0ReactionField(ISpace space, P2ReactionFieldDipole p) {
            super(0,space);
            this.potential = p;
            gradient = new IVectorMutable[0];
        }
        
        public double energy(IMoleculeList atoms) {
            double epsilon = potential.getDielectric();
            double cutoff = potential.getRange();
            DipoleSource dipoleSource = potential.getDipoleSource();
            double fac = 2*(epsilon-1)/(2*epsilon+1)/(cutoff*cutoff*cutoff);
            double u = 0;
            if (targetAtom != null) {
                IVector iDipole = dipoleSource.getDipole(targetAtom);
                u = -0.5 * fac * iDipole.squared();
            }
            else {
                IMoleculeList moleculeList = box.getMoleculeList();
                for (int i=0; i<moleculeList.getMoleculeCount(); i++) {
                    IVector iDipole = dipoleSource.getDipole(moleculeList.getMolecule(i));
                    u += -0.5 * fac * iDipole.squared();
                }
            }
            return u;
        }
        
        public void setBox(IBox newBox) {
            box = newBox;
        }
        
        public void setTargetMolecule(IMolecule atom) {
            if (atom == null) {
                targetAtom = null;
                return;
            }
            targetAtom = atom;
        }
        
        public void setTargetAtom(IAtom targetAtom) {
            throw new RuntimeException("Can't provide correction for an individual atom");
        }

        public double getRange() {
            return 0;
        }

        public IVector[] gradient(IMoleculeList atoms) {
            return gradient;
        }
        
        public IVector[] gradient(IMoleculeList atoms, Tensor pressureTensor) {
            return gradient(atoms);
        }
        
        public double virial(IMoleculeList atoms) {
            return 0;
        }

        protected final P2ReactionFieldDipole potential;
        protected final IVectorMutable[] gradient;
        protected IMolecule targetAtom;
        protected IBox box;

    }
}
