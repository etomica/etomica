/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.liquidLJ;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.normalmode.CoordinateDefinition;
import etomica.potential.Potential2;
import etomica.potential.Potential2Soft;
import etomica.potential.PotentialSoft;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;

/**
 * Class for spherical pair potential that computes lattice sums for various
 * truncations.  This class only pretends to implement the PotentialSoft
 * interface.  Only the energyVirialCut method should be called.  It will
 * return an instance of ReturnValue, which contains sums of various
 * quantities for all cutoffs passed to the constructor.
 *
 * This will compute self-contributions (interaction between an Atom and its
 * own image), but only if invoked with an atom pair holding two intsances of
 * that Atom (PotentialMaster classes don't know how to do this).
 *
 * @author Andrew Schultz
 */
 
public class Potential2SoftSphericalLSMultiLat extends Potential2 implements PotentialSoft {
   
    public Potential2SoftSphericalLSMultiLat(Space space, double[] rCut, Potential2Soft p2Soft, CoordinateDefinition coordinateDefinition) {
        super(space);
        gradient = new Vector[2];
        gradient[0] = space.makeVector();
        gradient[1] = space.makeVector();
        dr = space.makeVector();
        this.rCut2 = new double[rCut.length];
        for (int i=0; i<rCut2.length; i++) {
            rCut2[i] = rCut[i]*rCut[i];
        }
        this.p2Soft = p2Soft;
    	Lxyz = space.makeVector();
		drTmp = space.makeVector();
		drLat = space.makeVector();
		drLatTmp = space.makeVector();
		drA = space.makeVector();
		this.rCutMax = rCut[rCut.length-1];
		nShells = new int[space.D()];
		a0 = new double[space.D()];
		this.coordinateDefinition = coordinateDefinition;
		pTmp1 = space.makeVector();
		rv = new ReturnValue(rCut.length, space);
	}

    public double energy(IAtomList atoms) {
        return 0;
    }
    
    public ReturnValue energyVirialCut(IAtomList atoms) {
    	boolean isSelf = (atoms.get(1) == atoms.get(0));
        dr.Ev1Mv2(atoms.get(1).getPosition(),atoms.get(0).getPosition());
        drLat.E(coordinateDefinition.getLatticePosition(atoms.get(1)));
        drLat.ME(coordinateDefinition.getLatticePosition(atoms.get(0)));
        drTmp.E(drLat);
        boundary.nearestImage(drLat);
        dr.PE(drLat);
        dr.ME(drTmp);
        
        drTmp.Ev1Mv2(atoms.get(1).getPosition(), coordinateDefinition.getLatticePosition(atoms.get(1)));
        drA.E(drTmp);
        drTmp.Ev1Mv2(atoms.get(0).getPosition(), coordinateDefinition.getLatticePosition(atoms.get(0)));
        drA.ME(drTmp);
        
        for (int i=0; i<rCut2.length; i++) {
            rv.dadbSum[i] = rv.energySum[i] = rv.sum1[i] = rv.virialSum[i] = rv.pzxySum[i] = 0;
        }
        for(int nx = -nShells[0]; nx <= nShells[0]; nx++) {
        	Lxyz.setX(0, nx*a0[0]);
            for(int ny = -nShells[1]; ny <= nShells[1]; ny++) {
            	Lxyz.setX(1, ny*a0[1]);
                for(int nz = -nShells[2]; nz <= nShells[2]; nz++) {
                	Lxyz.setX(2, nz*a0[2]);
					drLatTmp.Ev1Pv2(drLat, Lxyz);
					double dr2Lat = drLatTmp.squared();
					if(dr2Lat > rCutMax*rCutMax ) continue;
                	boolean centerImage = (nx*nx+ny*ny+nz*nz == 0);
                	if(isSelf && centerImage) continue;
                    drTmp.Ev1Pv2(dr, Lxyz);
                    double dr2 = drTmp.squared();
                	double pu = p2Soft.u(dr2);
                	double dpu = p2Soft.du(dr2);
                    double p1 = -dpu/dr2*drTmp.dot(drLatTmp);
                    double dadb = -dpu/dr2*drTmp.dot(drA);
                    pTmp1.Ea1Tv1(-dpu / dr2, drTmp);
                    double pzxy = pTmp1.getX(2) * drTmp.getX(2) - (pTmp1.getX(0) * drTmp.getX(0) + pTmp1.getX(1) * drTmp.getX(1)) / 2;
                	if (isSelf) {
                	    pu *= 0.5;
                	    dpu *= 0.5;
                	    p1 *= 0.5;
                	    dadb *= 0.5;
                        pzxy *= 0.5;
                	}
                    for (int i=rCut2.length-1; i>=0; i--) {
                        if (dr2Lat > rCut2[i]) break;
                        rv.energySum[i] += pu;
                        rv.virialSum[i] += dpu;
                        rv.sum1[i] += p1;
                        rv.dadbSum[i] += dadb;
                        rv.pzxySum[i] += pzxy;
                    }
                }
            }
        }
        return rv;
    }

    /**
     * Virial of the pair as given by the du(double) method
     */
    public double virial(IAtomList atoms) {
        return 0;
    }
    
    
    /**
     * Gradient of the pair potential as given by the du(double) method.
     */
    public Vector[] gradient(IAtomList atoms) {
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

        
        for (int i=0; i<3; i++) {
            a0[i] = boundary.getBoxSize().getX(i);
            nShells[i] = (int) Math.ceil(rCutMax/a0[0] - 0.49999);
        }

    }

    protected final Vector[] gradient;
    protected Boundary boundary;
    protected final int[] nShells;
    protected final double[] a0;
    protected final Potential2Soft p2Soft;
    protected final Vector Lxyz;
    protected final Vector dr, drTmp, drLat, drA, drLatTmp;
    protected final Vector pTmp1;
    protected final double[] rCut2;
    protected final double rCutMax;
    protected final CoordinateDefinition coordinateDefinition;
    protected final ReturnValue rv;

    public class ReturnValue {
        public double[] energySum, virialSum, sum1, dadbSum, pzxySum;
        public ReturnValue(int n, Space space) {
            energySum = new double[n];
            virialSum = new double[n];
            sum1 = new double[n];
            dadbSum = new double[n];
            pzxySum = new double[n];
        }
    }
}//end of Potential2SoftSpherical
