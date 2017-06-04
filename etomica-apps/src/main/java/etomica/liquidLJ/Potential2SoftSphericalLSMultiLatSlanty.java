/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.liquidLJ;

import etomica.atom.IAtomList;
import etomica.space.Vector;
import etomica.box.Box;
import etomica.normalmode.CoordinateDefinition;
import etomica.potential.Potential2Soft;
import etomica.space.Space;

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
        
        for (int i=0; i<rCut2.length; i++) {
            rv.dadbSum[i] = rv.energySum[i] = rv.sum1[i] = rv.virialSum[i] = 0;
            rv.pSumXYZ1[i].E(0);
            rv.pSumXYZ2[i].E(0);
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
                    pTmp2.Ea1Tv1(-dpu/dr2, drTmp);
                    pTmp2.TE(drA);
                	if (isSelf) {
                	    pu *= 0.5;
                	    dpu *= 0.5;
                	    p1 *= 0.5;
                	    dadb *= 0.5;
                	    pTmp1.TE(0.5);
                        pTmp2.TE(0.5);
                	}
                    for (int i=rCut2.length-1; i>=0; i--) {
                        if (dr2Lat > rCut2[i]) break;
                        rv.energySum[i] += pu;
                        rv.virialSum[i] += dpu;
                        rv.sum1[i] += p1;
                        rv.dadbSum[i] += dadb;
                        rv.pSumXYZ1[i].PE(pTmp1);
                        rv.pSumXYZ2[i].PE(pTmp2);
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
