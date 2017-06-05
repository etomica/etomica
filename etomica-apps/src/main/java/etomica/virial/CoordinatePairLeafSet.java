/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.IAtomList;
import etomica.space.Vector;
import etomica.space.Space;
import etomica.util.Debug;

/**
 * @author David Kofke
 *
 * Class that holds a set of atom pairs.  Takes a list of atoms in its
 * constructor, and forms an instance of CoordinatePair for each pair formed in the
 * list.  Each CoordinatePair instance can be accessed via the getPair method.  The
 * CoordinatePairLeafSet has an ID, which changes whenever the atom positions change 
 * (after which reset() should be called).  Cluster.value() depends on a one-to-one
 * correspondence between this ID and the positions of atoms in the CoordinatePairLeafSet. 
 */

public class CoordinatePairLeafSet implements CoordinatePairSet {

    /**
     * Constructor for CoordinatePairLeafSet.
     * @param list The list of atoms for which the set of pairs is formed.
     */
    public CoordinatePairLeafSet(IAtomList list, Space space) {
        positions = new Vector[list.getAtomCount()];
        numAtoms = list.getAtomCount();
        r2 = new double[numAtoms*numAtoms];
        setAtoms(list);
        dr = space.makeVector();
    }
    
    /**
     * Returns atom pair for ith and jth atoms in set.
     */
    public double getr2(int i, int j) {
        if(Debug.ON && !(i<j)) throw new IllegalArgumentException("Error: i must be less than j");
        return r2[i*numAtoms+j];
    }

    private void setAtoms(IAtomList list) {
        for (int i=0; i<list.getAtomCount(); i++) {
            positions[i] = list.getAtom(i).getPosition();
        }
    }
    
    public void reset(long cPairID) {
        for(int i=0; i<numAtoms-1; i++) {
            Vector pos1 = positions[i];
            for(int j=i+1; j<numAtoms; j++) {
                dr.Ev1Mv2(pos1,positions[j]);
                r2[i*numAtoms+j] = dr.squared();
            }
        }
        ID = cPairID;
    }
    
    public void E(CoordinatePairLeafSet c) {
        for(int i=0; i<numAtoms-1; i++) {
            for(int j=i+1; j<numAtoms; j++) {
                r2[i*numAtoms+j] =  c.r2[i*numAtoms+j];
            }
        }
    }
    
    public long getID() {
        return ID;
    }
    
    protected final double[] r2;
    protected final Vector[] positions;
    protected final int numAtoms;
    protected final Vector dr;
    protected long ID;
}
