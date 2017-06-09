/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.IMolecule;
import etomica.atom.IMoleculeList;
import etomica.atom.IMoleculePositionDefinition;
import etomica.space.Vector;
import etomica.atom.MoleculePositionGeometricCenter;
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
public class CoordinatePairMoleculeSet implements CoordinatePairSet {

    /**
     * Constructor for CoordinatePairLeafSet.
     * @param list The list of atoms for which the set of pairs is formed.
     */
    public CoordinatePairMoleculeSet(IMoleculeList list, Space space) {
        numAtoms = list.getMoleculeCount();
        atoms = new IMolecule[numAtoms];
        r2 = new double[numAtoms*numAtoms];
        setAtoms(list);
        dr = space.makeVector();
        iPosition = space.makeVector();
        positionDefinition = new MoleculePositionGeometricCenter(space);
    }
    
    
    public IMoleculePositionDefinition getPositionDefinition() {
		return positionDefinition;
	}


	public void setPositionDefinition(IMoleculePositionDefinition positionDefinition) {
		this.positionDefinition = positionDefinition;
	}


	/**
     * Returns atom pair for ith and jth atoms in set.
     */
    public double getr2(int i, int j) {
        if(Debug.ON && !(i<j)) throw new IllegalArgumentException("Error: i must be less than j");
        return r2[i*numAtoms+j];
    }

    private void setAtoms(IMoleculeList list) {
        for (int i=0; i<list.getMoleculeCount(); i++) {
            atoms[i] = list.getMolecule(i);
        }
    }
    
    public void reset(long cPairID) {
        for(int i=0; i<numAtoms-1; i++) {
            IMolecule iAtom = atoms[i];
            iPosition.E(positionDefinition.position(iAtom));
            for(int j=i+1; j<numAtoms; j++) {
                IMolecule jAtom = atoms[j];
                Vector jPosition = positionDefinition.position(jAtom);
                dr.Ev1Mv2(iPosition, jPosition);
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
    protected final IMolecule[] atoms;
    protected final int numAtoms;
    protected final Vector dr;
    protected final Vector iPosition;
    protected long ID;
    protected IMoleculePositionDefinition positionDefinition;
}
