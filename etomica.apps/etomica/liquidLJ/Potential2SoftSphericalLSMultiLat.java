/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.liquidLJ;

import etomica.api.IAtomList;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.normalmode.CoordinateDefinition;
import etomica.potential.Potential2;
import etomica.potential.Potential2Soft;
import etomica.potential.PotentialSoft;
import etomica.space.ISpace;
import etomica.space.Tensor;

/**
 * Methods for a soft (non-impulsive), spherically-symmetric pair potential.
 * Subclasses must provide concrete definitions for the energy (method
 * u(double)) and its derivatives.
 * 
 * @author David Kofke
 */
 
public class Potential2SoftSphericalLSMultiLat extends Potential2 implements PotentialSoft {
   
    public Potential2SoftSphericalLSMultiLat(ISpace space, double[] rCut, Potential2Soft p2Soft, CoordinateDefinition coordinateDefinition) {
        super(space);
        gradient = new IVectorMutable[2];
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
		sums = new double[4][rCut.length];
		this.coordinateDefinition = coordinateDefinition;
	}

    public double energy(IAtomList atoms) {
        return 0;
    }
    
    public double[][] energyVirialCut(IAtomList atoms) {
    	boolean isSelf = (atoms.getAtom(1) == atoms.getAtom(0));
        dr.Ev1Mv2(atoms.getAtom(1).getPosition(),atoms.getAtom(0).getPosition());
        drLat.E(coordinateDefinition.getLatticePosition(atoms.getAtom(1)));
        drLat.ME(coordinateDefinition.getLatticePosition(atoms.getAtom(0)));
        drTmp.E(drLat);
        boundary.nearestImage(drLat);
        dr.PE(drLat);
        dr.ME(drTmp);
        
        drTmp.Ev1Mv2(atoms.getAtom(1).getPosition(), coordinateDefinition.getLatticePosition(atoms.getAtom(1)));
        drA.E(drTmp);
        drTmp.Ev1Mv2(atoms.getAtom(0).getPosition(), coordinateDefinition.getLatticePosition(atoms.getAtom(0)));
        drA.ME(drTmp);
        
        for (int i=rCut2.length-1; i>=0; i--) {
            for (int j=0; j<4; j++) {
                sums[j][i] = 0;
            }
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
                	if (isSelf) {
                	    pu *= 0.5;
                	    dpu *= 0.5;
                	    p1 *= 0.5;
                	    dadb *= 0.5;
                	}
                    for (int i=rCut2.length-1; i>=0; i--) {
                        if (dr2Lat > rCut2[i]) break;
                        sums[0][i] += pu;
                        sums[1][i] += dpu;
                        sums[2][i] += p1;
                        sums[3][i] += dadb;
                    }
                }
            }
        }
        return sums;
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
    public IVector[] gradient(IAtomList atoms) {
        return gradient;
    }
    
    public IVector[] gradient(IAtomList atoms, Tensor pressureTensor) {
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
    
    public void setBox(IBox box) {
        boundary = box.getBoundary();
        p2Soft.setBox(box);

        
        for (int i=0; i<3; i++) {
            a0[i] = boundary.getBoxSize().getX(i);
            nShells[i] = (int) Math.ceil(rCutMax/a0[0] - 0.49999);
        }

    }

    protected final IVectorMutable[] gradient;
    protected IBoundary boundary;
    protected final int[] nShells;
    protected final double[] a0;
    protected final Potential2Soft p2Soft;
    protected final IVectorMutable Lxyz;
    protected final IVectorMutable dr, drTmp, drLat, drA, drLatTmp;
    protected final double[] rCut2;
    protected final double rCutMax;
    protected final double[][] sums;
    protected final CoordinateDefinition coordinateDefinition;

}//end of Potential2SoftSpherical
