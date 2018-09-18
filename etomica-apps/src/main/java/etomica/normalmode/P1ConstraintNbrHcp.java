/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.AtomArrayList;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.potential.IPotentialAtomic;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;

public class P1ConstraintNbrHcp implements IPotentialAtomic {

    /**
     * 
     */
    private static final long serialVersionUID = 1L;

    // this could take a NeighborListManager to try to speed up finding neighbors
    public P1ConstraintNbrHcp(Space space, double neighborDistance, Box box) {
        boundary = box.getBoundary();

        neighborRadiusSq = neighborDistance*neighborDistance;

        IAtomList list = box.getLeafList();

        //Check for neighboring sites
        drj = space.makeVector();
        drk = space.makeVector();
        neighborAtoms = new int[list.size()][12];
        AtomArrayList tmpList = new AtomArrayList(12);

        for (int i = 0; i<list.size(); i++) {
            IAtom atomi = list.get(i);
            tmpList.clear();
            for (int j = 0; j<list.size(); j++) {
                if (i==j) continue;
                IAtom atomj = list.get(j);
                drj.Ev1Mv2(atomi.getPosition(), atomj.getPosition());
                boundary.nearestImage(drj);
                if (drj.squared() < neighborRadiusSq*1.01) {
                    tmpList.add(atomj);
                }
            }
            for (int j=0; j<12; j++) {
                neighborAtoms[i][j] = tmpList.get(j).getLeafIndex();
            }
        }
    }

    public int nBody() {
        return 1;
    }

    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public void setBox(Box box) {
        leafList = box.getLeafList();
    }


    /**
     * Returns sum of energy for all triplets containing the given atom
     */
	public double energy(IAtomList atoms) {
	    IAtom atom = atoms.get(0);
	    double u = energyi(atom);
	    if (u == Double.POSITIVE_INFINITY) {
	        return u;
	    }

	    return u;
	}

	/**
	 * Returns the energy for atom i due to requirements for i and its
	 * neighbors (but not for i as a neighbor)
	 */
	public double energyi(IAtom atom) {

	    Vector posAtom = atom.getPosition();

	    int atomIndex = atom.getLeafIndex();
	    int[] list = neighborAtoms[atomIndex];
	    for (int i=0; i<12; i++) {
	        IAtom atomj = leafList.get(list[i]);
	        drj.Ev1Mv2(posAtom, atomj.getPosition());
	        boundary.nearestImage(drj);
	        if (drj.squared() > neighborRadiusSq*3.0) {
	        	p1Counter++;
	            return Double.POSITIVE_INFINITY;
	        }
	   
	    }
	    return 0;
	}
	
	public int getp1Counter(){
		return p1Counter;
	}
	

	protected final int[][] neighborAtoms;
	protected final Vector drj, drk;
	protected double neighborRadiusSq;
	protected final Boundary boundary;
	protected IAtomList leafList;
	protected int p1Counter=0;
}
