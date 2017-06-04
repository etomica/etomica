/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtom;
import etomica.api.IBoundary;
import etomica.box.Box;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.IPotentialMolecular;
import etomica.space.Vector;
import etomica.atom.DipoleSource;
import etomica.atom.IAtomPositionDefinition;
import etomica.space.Space;
import etomica.space.Tensor;

public class P2ReactionFieldDipole extends PotentialMolecular implements PotentialMolecularSoft, IPotentialMolecularSecondDerivative {

    public P2ReactionFieldDipole(Space space, IAtomPositionDefinition positionDefinition) {
        super(2, space);
        this.positionDefinition = positionDefinition;
        iDipole = space.makeVector();
        cavityDipole = space.makeVector();
        dr = space.makeVector();
        dr1 = space.makeVector();
        gradientAndTorque = new Vector[2][2];
        gradientAndTorque[0][0] = space.makeVector();
        gradientAndTorque[0][1] = space.makeVector();
        gradientAndTorque[1][0] = space.makeVector();
        gradientAndTorque[1][1] = space.makeVector();
        secondDerivative = new Tensor [3];
        secondDerivative[0] = space.makeTensor();
		secondDerivative[1] = space.makeTensor();
		secondDerivative[2] = space.makeTensor();
		a = new Vector[3];
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
        if (epsilon < Double.POSITIVE_INFINITY) {
            fac = 2*(epsilon-1)/(2*epsilon+1)/(cutoff2*cutoff);
        }
        else {
            fac = 1 / (cutoff2*cutoff);
        }
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
            if (epsilon < Double.POSITIVE_INFINITY) {
                fac = 2*(epsilon-1)/(2*epsilon+1)/(cutoff2*cutoff);
            }
            else {
                fac = 1 / (cutoff2*cutoff);
            }
        }
    }

    public void setBox(Box box) {
        boundary = box.getBoundary();
    }

    public double energy(IMoleculeList molecules) {
    	IMolecule molecule0 = molecules.getMolecule(0);
		IMolecule molecule1 = molecules.getMolecule(1);
		dr.E(positionDefinition.position(molecule1));
		dr.ME(positionDefinition.position(molecule0));
		boundary.nearestImage(dr);//rij
		double r2 = dr.squared();
		if (r2 > cutoff*cutoff) return 0;
		
        iDipole.E(dipoleSource.getDipole(molecules.getMolecule(0)));
        double idotj = iDipole.dot(dipoleSource.getDipole(molecules.getMolecule(1)));
        return -fac*idotj;
    }

    public Vector[][] gradientAndTorque(IMoleculeList molecules) {
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
		iDipole.E(dipoleSource.getDipole(molecule0));
		Vector jDipole = space.makeVector();
		jDipole.E(dipoleSource.getDipole(molecule1));
		
		double exi = iDipole.getX(0);//ei and ej is the dipole orientation with mu
		double eyi = iDipole.getX(1);
		double ezi = iDipole.getX(2);
		double exj = jDipole.getX(0);
		double eyj = jDipole.getX(1);
		double ezj = jDipole.getX(2);
		
		Vector deidxi = space.makeVector();
		Vector deidyi = space.makeVector();
		Vector deidzi = space.makeVector();
		Vector dejdxj = space.makeVector();
		Vector dejdyj = space.makeVector();
		Vector dejdzj = space.makeVector();
		
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
		
		a[0].E(iDipole);
		a[0].XE(dejdxj);
		a[1].E(iDipole);
		a[1].XE(dejdyj);
		a[2].E(iDipole);
		a[2].XE(dejdzj);
		secondDerivative[0].E(a);
		
		secondDerivative[0].TE(-fac);//ij and ji are the same
		
		a[0].E(deidxi);
		a[0].XE(jDipole);
		a[1].E(deidyi);
		a[1].XE(jDipole);
		a[2].E(deidzi);
		a[2].XE(jDipole);
		secondDerivative[1].E(a);
		
		secondDerivative[1].TE(-fac);//ii
		
		a[0].E(dejdxj);
		a[0].XE(iDipole);
		a[1].E(dejdyj);
		a[1].XE(iDipole);
		a[2].E(dejdzj);
		a[2].XE(iDipole);
		secondDerivative[2].E(a);
		
		
		secondDerivative[2].TE(-fac);//jj 
    	return secondDerivative; 
    }
    
    

    public Vector[] gradient(IMoleculeList atoms) {
        return gradientAndTorque[0];
    }

    public Vector[] gradient(IMoleculeList atoms, Tensor pressureTensor) {
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
    protected final Vector iDipole, cavityDipole;
    protected final Vector dr ,dr1;
    protected DipoleSource dipoleSource;
    protected IBoundary boundary;
    protected double cutoff2, cutoff;
    protected double epsilon;
    protected final Vector[][] gradientAndTorque;
    protected final Tensor[] secondDerivative;
    protected final Vector[] a;
    protected double fac;
    
    /**
     * A 0-body potential that should be added along with this potential.  The
     * 0-body potential includes the effective self-interaction of the
     * molecules (the molecule induces a dipole in the surrounding fluid, which
     * has an interaction energy with the molecule).  This part of the
     * potential does not result in a gradient or torque on the molecule and is
     * independent of position or orientation.
     */
    public static class P0ReactionField extends PotentialMolecular implements IPotential0Lrc, PotentialMolecularSoft {

        public P0ReactionField(Space space, P2ReactionFieldDipole p) {
            super(0,space);
            this.potential = p;
            gradient = new Vector[0];
        }
        
        public double energy(IMoleculeList atoms) {
            double epsilon = potential.getDielectric();
            double cutoff = potential.getRange();
            DipoleSource dipoleSource = potential.getDipoleSource();
            double fac = 2*(epsilon-1)/(2*epsilon+1)/(cutoff*cutoff*cutoff);
            double u = 0;
            if (targetAtom != null) {
                Vector iDipole = dipoleSource.getDipole(targetAtom);
                u = -0.5 * fac * iDipole.squared();
            }
            else {
                IMoleculeList moleculeList = box.getMoleculeList();
                for (int i=0; i<moleculeList.getMoleculeCount(); i++) {
                    Vector iDipole = dipoleSource.getDipole(moleculeList.getMolecule(i));
                    u += -0.5 * fac * iDipole.squared();
                }
            }
            return u;
        }
        
        public void setBox(Box newBox) {
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

        public Vector[] gradient(IMoleculeList atoms) {
            return gradient;
        }
        
        public Vector[] gradient(IMoleculeList atoms, Tensor pressureTensor) {
            return gradient(atoms);
        }
        
        public double virial(IMoleculeList atoms) {
            return 0;
        }

        protected final P2ReactionFieldDipole potential;
        protected final Vector[] gradient;
        protected IMolecule targetAtom;
        protected Box box;

    }
}
