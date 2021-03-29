/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;

/**
 * Hard sphere molecule with a dipole sitting at the center.
 *
 * @author Weisong & Shu
 * date:May 2015
 * 
 */
public class P2HSDipole extends PotentialMolecular implements IPotentialMolecularSecondDerivative {

	public P2HSDipole(Space space, double sigma, double dipole) {
		this(space, 1,dipole,Double.POSITIVE_INFINITY);
	}

	public P2HSDipole(Space space, double sigma, double dipole, double rCut) {
		super(2, space);
		setSigma(sigma);
		gradient = new Vector[2];
		gradient[0] = space.makeVector();
		gradient[1] = space.makeVector();
		torque = new Vector[2];
		torque[0] = space.makeVector();
		torque[1] = space.makeVector();
		secondDerivative = new Tensor[3];
		this.secondDerivative[0] = space.makeTensor();
		this.secondDerivative[1] = space.makeTensor();
		this.secondDerivative[2] = space.makeTensor();
		a = new Vector[3];
		a[0] = space.makeVector();
		a[1] = space.makeVector();
		a[2] = space.makeVector();
		
		dr = space.makeVector();
		drunit = space.makeVector();
		runit = space.makeVector();
		work = space.makeVector();
		gradientAndTorque = new Vector[][]{gradient,torque};

		setDipole(dipole);
		this.rCut = rCut;
	}


	public void setBox(Box box) {
		boundary = box.getBoundary();
	}

	public double getRange() {
		return Double.POSITIVE_INFINITY;
	}

	public double energy(IMoleculeList pair){
		IMolecule molecule1 = pair.get(0);
		IMolecule molecule2 = pair.get(1);
		IAtomList atomList1 = molecule1.getChildList();
		IAtomList atomList2 = molecule2.getChildList();

		IAtomOriented atom1 = (IAtomOriented)atomList1.get(0);
		IAtomOriented atom2 = (IAtomOriented)atomList2.get(0);
		// dr is a r1-r2
		dr.Ev1Mv2(atom1.getPosition(), atom2.getPosition());
		boundary.nearestImage(dr);
		double r2 = dr.squared();

		if(r2 > rCut*rCut){
			return 0;
		}

		if(r2 < sigma2) { // hard core
			return Double.POSITIVE_INFINITY;
		}
		// normalize dr, the vector between the molecules
		dr.normalize();

		// v1 (unit vector) is the orientation of molecule 1: dipole1 direction
		Vector v1 = atom1.getOrientation().getDirection();
		// v2 (unit vector) is the orientation of molecule 2: dipole1 direction
		Vector v2 = atom2.getOrientation().getDirection();

		// cos(dipole 1 and dipole 2)=cos(v1 and v2)
		double cos_D1_D2 = v1.dot(v2);
		//cos(dipole 1 and r12)
		double cos_D1_r = v1.dot(dr);
		//cos(r12 and dipole 2)
		double cos_r_D2=dr.dot(v2);
		double r12Magnitude = Math.sqrt(r2);
		double ener = dipole * dipole *  (cos_D1_D2 - 3.0  * cos_D1_r* cos_r_D2);
		ener = ener/r12Magnitude /r2 ;
		return ener;
	}

	public double getSigma() {return sigma;}

	public final void setSigma(double s) {
		sigma = s;
		sigma2 = s * s;
	}

	public void setDipole(double moment){
		dipole=moment;
	}
	
    public Vector[] gradient(IMoleculeList pair, Tensor pressureTensor) {
        return gradient(pair);
    }

    public Vector[] gradient(IMoleculeList pair) {
        // do extra work to calculate torque
        gradientAndTorque(pair);
        return gradient;
    }
    
    public double virial(IMoleculeList atoms) {
        gradient(atoms);
    	IMolecule molecule1 = atoms.get(0);
		IMolecule molecule2 = atoms.get(1);
		IAtomList atomList1 = molecule1.getChildList();
		IAtomList atomList2 = molecule2.getChildList();

		IAtomOriented atom1 = (IAtomOriented)atomList1.get(0);
		IAtomOriented atom2 = (IAtomOriented)atomList2.get(0);

        // LJ contributation

        dr.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
        boundary.nearestImage(dr);
        double v = gradient[1].dot(dr);
        if (Double.isInfinite(v)) {
            throw new RuntimeException("oops "+v);
        }
        return gradient[1].dot(dr);
    }

	public Vector[][] gradientAndTorque(IMoleculeList atoms) {
    	IMolecule molecule1 = atoms.get(0);
		IMolecule molecule2 = atoms.get(1);
		IAtomList atomList1 = molecule1.getChildList();
		IAtomList atomList2 = molecule2.getChildList();

		IAtomOriented atom1 = (IAtomOriented)atomList1.get(0);
		IAtomOriented atom2 = (IAtomOriented)atomList2.get(0);

		dr.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
		boundary.nearestImage(dr);
		double r2 = dr.squared();
		if (r2 > rCut*rCut) {
			gradient[0].E(0);
			gradient[1].E(0);
			torque[0].E(0);
			torque[1].E(0);
			return gradientAndTorque;
		}

		gradient[0].E(0);
		double momentSq = dipole*dipole;
		if(momentSq!=0.0){
			double s2 = 1/r2;
			double s1 = Math.sqrt(s2);
			// normalize dr, the vector between the molecules
			drunit.E(dr);
			drunit.TE(s1);					

			// v1 is the orientation of molecule 1
			Vector v1 = atom1.getOrientation().getDirection();

			// v2 is the orientation of molecule 2
			Vector v2 = atom2.getOrientation().getDirection();

			double fac = momentSq * (s2*s1 );
			double dfac = momentSq * (3*s2*s2*s1 );
			double udd = v1.dot(v2) - 3.0*v1.dot(drunit)*v2.dot(drunit);
			work.Ea1Tv1(v2.dot(drunit), v1);
			work.PEa1Tv1(v1.dot(drunit), v2);
			work.PEa1Tv1(-2*v1.dot(drunit)*v2.dot(drunit)*s1, dr);
			work.TE(3.0*s1*fac);
			gradient[0].PE(work);
			gradient[0].PEa1Tv1(dfac*udd, dr);

			work.E(v1);
			work.XE(v2);
			torque[0].E(v1);
			torque[0].XE(drunit);
			torque[0].TE(3.0*v2.dot(drunit));
			torque[0].ME(work);
			torque[0].TE(fac);

			torque[1].E(v2);
			torque[1].XE(drunit);
			torque[1].TE(3.0*v1.dot(drunit));
			torque[1].PE(work);
			torque[1].TE(fac);

			
//			System.out.println("dr = " + dr);
//			System.out.println("vi = " + v1);
//			System.out.println("vj = " + v2);
//			System.out.println("torque[0] = " + torque[0]);   //TODO debug only
//			System.out.println("torque[1] = " + torque[1]);   //TODO debug only
		}

		// pairwise additive, so
		gradient[1].Ea1Tv1(-1,gradient[0]);

		return gradientAndTorque;
	}
	
	public Tensor[] secondDerivative(IMoleculeList molecules){
		IMolecule molecule0 = molecules.get(0);
		IMolecule molecule1 = molecules.get(1);
		IAtomList atomList0 = molecule0.getChildList();
		IAtomList atomList1 = molecule1.getChildList();
		IAtomOriented atom0 = (IAtomOriented)atomList0.get(0);
		IAtomOriented atom1 = (IAtomOriented)atomList1.get(0);
		Vector pos0 = atom1.getPosition();
		Vector pos1 = atom0.getPosition();
		Vector ei =  atom0.getOrientation().getDirection();
		Vector ej =  atom1.getOrientation().getDirection();
		
		double exi = ei.getX(0);//ei and ej is the dipole orientation
		double eyi = ei.getX(1);
		double ezi = ei.getX(2);
		double exj = ej.getX(0);
		double eyj = ej.getX(1);
		double ezj = ej.getX(2);
		
		Vector deidxi = space.makeVector();
		Vector deidyi = space.makeVector();
		Vector deidzi = space.makeVector();
		Vector dejdxj = space.makeVector();
		Vector dejdyj = space.makeVector();
		Vector dejdzj = space.makeVector();
		
		double [] dejdxjD = {0,-ezj,eyj};
		double [] dejdyjD = {ezj,0,-exj};
		double [] dejdzjD = {-eyj,exj,0};
		dejdxj.E(dejdxjD);
		dejdyj.E(dejdyjD);
		dejdzj.E(dejdzjD);
		
		dr.Ev1Mv2(pos1, pos0);
		boundary.nearestImage(dr);//rij
		
		
//		System.out.println("rij_P2 = " + dr);// debug only
		
		
		double r2 = dr.squared();
		secondDerivative[0].E(0);
		secondDerivative[1].E(0);
		secondDerivative[2].E(0);
		if (r2 > rCut*rCut) return secondDerivative;
		double r = Math.sqrt(r2);
		runit.Ea1Tv1(1.0/r, dr);
		double coeff = dipole*dipole/r2/r;
		
		//ei cross dejdxj - 3.0*(dejdxj.runit)*(ei cross runite)
		dr.E(ei);
		dr.XE(dejdxj);
		a[0].E(ei);
		a[0].XE(runit);
		a[0].TE(-3.0*dejdxj.dot(runit));
		a[0].PE(dr);
		
		dr.E(ei);
		dr.XE(dejdyj);
		a[1].E(ei);
		a[1].XE(runit);
		a[1].TE(-3.0*dejdyj.dot(runit));
		a[1].PE(dr);
		
		dr.E(ei);
		dr.XE(dejdzj);
		a[2].E(ei);
		a[2].XE(runit);
		a[2].TE(-3.0*dejdzj.dot(runit));
		a[2].PE(dr);
		
		secondDerivative[0].E(a);
		secondDerivative[0].TE(coeff);//ij
		
//		System.out.println("ei = " + ei);// debug only
//		System.out.println("ej = " + ej);// debug only
//		System.out.println("dudij = \n" + secondDerivative[0]);// debug only
		
		
		
		double [] deidxiD = {0,-ezi,eyi};
		double [] deidyiD = {ezi,0,-exi};
		double [] deidziD = {-eyi,exi,0};
		deidxi.E(deidxiD);
		deidyi.E(deidyiD);
		deidzi.E(deidziD);
		
//		deidxi cross ej - 3.0*(deidxi cross runite)*(ej.runit)
		dr.E(deidxi);
		dr.XE(ej);
		a[0].E(deidxi);
		a[0].XE(runit);
		a[0].TE(-3.0*ej.dot(runit));
		a[0].PE(dr);
		
		dr.E(deidyi);
		dr.XE(ej);
		a[1].E(deidyi);
		a[1].XE(runit);
		a[1].TE(-3.0*ej.dot(runit));
		a[1].PE(dr);
		
		dr.E(deidzi);
		dr.XE(ej);
		a[2].E(deidzi);
		a[2].XE(runit);
		a[2].TE(-3.0*ej.dot(runit));
		a[2].PE(dr);
		
		
		secondDerivative[1].E(a);
		secondDerivative[1].TE(coeff);//ii
		
//		System.out.println("dudii = \n" + secondDerivative[1]);//TODO debug only
		
		//dejdxj cross ei -3.0*(ei.runit)*(dejdxj cross runit)    
		dr.E(dejdxj);
		dr.XE(ei);
		a[0].E(dejdxj);
		a[0].XE(runit);
		a[0].TE(-3.0*ei.dot(runit));
		a[0].PE(dr);
		
		dr.E(dejdyj);
		dr.XE(ei);
		a[1].E(dejdyj);
		a[1].XE(runit);
		a[1].TE(-3.0*ei.dot(runit));
		a[1].PE(dr);
		
		dr.E(dejdzj);
		dr.XE(ei);
		a[2].E(dejdzj);
		a[2].XE(runit);
		a[2].TE(-3.0*ei.dot(runit));
		a[2].PE(dr);
		
		secondDerivative[2].E(a);
		secondDerivative[2].TE(coeff);//jj
		
		
//		System.out.println("dudjj = \n" + secondDerivative[2]);//TODO debug only
//		System.exit(2);// TODO bebug only
        return secondDerivative;
        
    }
	

	private static final long serialVersionUID = 1L;
	private double sigma , sigma2;
	private Boundary boundary;
	private final Vector dr,drunit,work,runit;
	private double dipole;
	private double rCut;
	private final Vector[] gradient;
	protected final Vector[] torque;
	protected final Vector[][] gradientAndTorque;
	protected final Vector[] a;
	protected final Tensor[] secondDerivative;		

}
