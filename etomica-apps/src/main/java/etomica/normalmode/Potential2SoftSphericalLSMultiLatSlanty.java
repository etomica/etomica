/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.potential.Potential2Soft;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Methods for a soft (non-impulsive), spherically-symmetric pair potential.
 * Subclasses must provide concrete definitions for the energy (method
 * u(double)) and its derivatives.
 * 
 * @author David Kofke
 */
 
public class Potential2SoftSphericalLSMultiLatSlanty extends Potential2SoftSphericalLSMultiLat {
   
    public Potential2SoftSphericalLSMultiLatSlanty(Space space, double[] rCut, Potential2Soft p2Soft, CoordinateDefinition coordinateDefinition, int[] nShells) {
        super(space, rCut, p2Soft, coordinateDefinition);
        xLxyz = space.makeVector();
        yLxyz = space.makeVector();
        zLxyz = space.makeVector();
		xLat = space.makeVector();
		yLat = space.makeVector();
		zLat = space.makeVector();
		System.arraycopy(nShells, 0, this.nShells, 0, 3);
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
        	xLxyz.Ea1Tv1(nx, xLat);
            for(int ny = -nShells[1]; ny <= nShells[1]; ny++) {
                yLxyz.E(xLxyz);
                yLxyz.PEa1Tv1(ny, yLat);
                for(int nz = -nShells[2]; nz <= nShells[2]; nz++) {
                    zLxyz.E(yLxyz);
                    zLxyz.PEa1Tv1(nz, zLat);
					drLatTmp.Ev1Pv2(drLat, zLxyz);
					double dr2Lat = drLatTmp.squared();
					if(dr2Lat > rCutMax*rCutMax ) continue;
                	boolean centerImage = (nx*nx+ny*ny+nz*nz == 0);
                	if(isSelf && centerImage) continue;
                    drTmp.Ev1Pv2(dr, zLxyz);
                    double dr2 = drTmp.squared();
                	double pu = p2Soft.u(dr2);
                	double dpu = p2Soft.du(dr2);
                    double p1 = -dpu/dr2*drTmp.dot(drLatTmp);
                    pTmp1.Ea1Tv1(-dpu/dr2, drTmp);
                    pTmp1.TE(drLatTmp);
                    double dadb = -dpu/dr2*drTmp.dot(drA);
                    pTmp1.Ea1Tv1(-dpu / dr2, drTmp);
                    double pzxy = pTmp1.getX(2) * dr.getX(2) - (pTmp1.getX(0) * dr.getX(0) + pTmp1.getX(1) * dr.getX(1)) / 2;
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
    
    public void setBox(Box box) {
        boundary = box.getBoundary();
        p2Soft.setBox(box);

        xLat.E(boundary.getEdgeVector(0));
        yLat.E(boundary.getEdgeVector(1));
        zLat.E(boundary.getEdgeVector(2));
    }

    protected final Vector xLxyz, yLxyz, zLxyz;
    protected final Vector xLat, yLat, zLat;
}//end of Potential2SoftSpherical
