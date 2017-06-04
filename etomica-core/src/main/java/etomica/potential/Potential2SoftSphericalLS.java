/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtomList;
import etomica.api.IBoundary;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.space.Space;
import etomica.space.Tensor;

/**
 * Methods for a soft (non-impulsive), spherically-symmetric pair potential.
 * Subclasses must provide concrete definitions for the energy (method
 * u(double)) and its derivatives.
 * 
 * @author David Kofke
 */
 
public class Potential2SoftSphericalLS extends Potential2 implements PotentialSoft{
   
    public Potential2SoftSphericalLS(Space space, double rCut, double[] a0, Potential2Soft p2Soft) {
         super(space);
        gradient = new Vector[2];
        gradient[0] = space.makeVector();
        gradient[1] = space.makeVector();
        dr = space.makeVector();
        this.rCut2 = rCut*rCut;
        this.a0 = a0;
        this.p2Soft = p2Soft;
    	Lxyz = space.makeVector();
		drtmp = space.makeVector();
		nShells = new int[] {(int) Math.ceil(rCut/a0[0] - 0.49999), (int) Math.ceil(rCut/a0[1] - 0.49999), (int) Math.ceil(rCut/a0[2] - 0.49999)};
	}
        
    public double energy(IAtomList atoms) {
    	boolean isSelf = (atoms.getAtom(1) == atoms.getAtom(0));
		double u_LJ = 0;
        dr.Ev1Mv2(atoms.getAtom(1).getPosition(),atoms.getAtom(0).getPosition());
        boundary.nearestImage(dr);
        for(int nx = -nShells[0]; nx <= nShells[0]; nx++) {
        	Lxyz.setX(0, nx*a0[0]);
            for(int ny = -nShells[1]; ny <= nShells[1]; ny++) {
            	Lxyz.setX(1, ny*a0[1]);
                for(int nz = -nShells[2]; nz <= nShells[2]; nz++) {
                	Lxyz.setX(2, nz*a0[2]);
					drtmp.Ev1Pv2(dr, Lxyz);
					double dr2 = drtmp.squared();
					if(dr2 > rCut2 ) continue;
                	boolean centerImage = (nx*nx+ny*ny+nz*nz == 0);
                	if(isSelf && centerImage) continue;
                	u_LJ += (isSelf ? 0.5 : 1.0)*p2Soft.u(dr2);
                }
            }
        }
        return u_LJ;
    }
    
    /**
     * Virial of the pair as given by the du(double) method
     */
    public double virial(IAtomList atoms) {
        double tmpVir = 0;
        dr.Ev1Mv2(atoms.getAtom(1).getPosition(),atoms.getAtom(0).getPosition());
        boundary.nearestImage(dr);
		for(int nx = -nShells[0]; nx <= nShells[0]; nx++) {
        	Lxyz.setX(0, nx*a0[0]);
            for(int ny = -nShells[1]; ny <= nShells[1]; ny++) {
            	Lxyz.setX(1, ny*a0[1]);
                for(int nz = -nShells[2]; nz <= nShells[2]; nz++) {
                	Lxyz.setX(2, nz*a0[2]);
					drtmp.Ev1Pv2(dr, Lxyz);
					tmpVir += p2Soft.du(drtmp.squared());
                }
            }
        }
      return tmpVir;
    }
    
    
    /**
     * Gradient of the pair potential as given by the du(double) method.
     */
    public Vector[] gradient(IAtomList atoms) {
    	boolean isSelf = (atoms.getAtom(1) == atoms.getAtom(0));
        dr.Ev1Mv2(atoms.getAtom(1).getPosition(),atoms.getAtom(0).getPosition());
        boundary.nearestImage(dr);
        gradient[0].E(0);
        gradient[1].E(0);
        for(int nx = -nShells[0]; nx <= nShells[0]; nx++) {
        	Lxyz.setX(0, nx*a0[0]);
            for(int ny = -nShells[1]; ny <= nShells[1]; ny++) {
            	Lxyz.setX(1, ny*a0[1]);
                for(int nz = -nShells[2]; nz <= nShells[2]; nz++) {
                	boolean nonCenterImage = (nx*nx+ny*ny+nz*nz > 0);
                	if(isSelf && nonCenterImage) continue;
                	Lxyz.setX(2, nz*a0[2]);
					drtmp.Ev1Pv2(dr, Lxyz);
					double dr2 = drtmp.squared();
					if(dr2 > rCut2 ) continue;
			        gradient[1].PEa1Tv1(p2Soft.du(drtmp.squared())/drtmp.squared(),drtmp);
                }
            }
        }
        gradient[0].PEa1Tv1(-1,gradient[1]);
        return gradient;
    }
    
    public Vector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        gradient(atoms);
        pressureTensor.PEv1v2(gradient[0],dr);
        return gradient;
    }
    
    
    /**
     * Returns infinity.  May be overridden to define a finite-ranged potential.
     */
    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public void setBox(Box box) {
        boundary = box.getBoundary();
        p2Soft.setBox(box);
    }

    protected final Vector[] gradient;
    protected IBoundary boundary;
    protected final int[] nShells;
    protected final double[] a0;
    protected final Potential2Soft p2Soft;
    protected final Vector Lxyz;
    protected final Vector dr;
    protected final Vector drtmp;
    protected final double rCut2;
    

}//end of Potential2SoftSpherical
