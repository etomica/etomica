/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.IMolecule;
import etomica.atom.IMoleculeList;
import etomica.atom.AtomPair;

/**
 * Iterator that returns intramolecular pairs of atoms with specific indices.
 * 
 * @author Andrew Schultz
 */
public class ApiIndexList implements AtomsetIteratorBasisDependent {

    /**
     * Constructs iterator without defining atoms in pair.
     */
    public ApiIndexList(int[][] index) {
        pair = new AtomPair();
        this.index = index;
        int max = -1;
        for (int i=0; i<index.length; i++) {
        	for (int j=0; j<2; j++) {
        		if (index[i][j] > max) max=index[i][j];
        	}
        }
        partners = new int[max+1][2][];
        int[][] nPartners = new int[partners.length][2];
        for (int i=0; i<index.length; i++) {
        	nPartners[index[i][0]][0]++;
        	nPartners[index[i][1]][1]++;
        }
        for (int i=0; i<partners.length; i++) {
            partners[i][0] = new int[nPartners[i][0]];
            partners[i][1] = new int[nPartners[i][1]];
        }
        nPartners = new int[partners.length][2];
        for (int i=0; i<index.length; i++) {
        	int idx0 = index[i][0];
        	int idx1 = index[i][1];
        	partners[idx0][0][nPartners[idx0][0]] = idx1;
            partners[idx1][1][nPartners[idx1][1]] = idx0;
            nPartners[idx0][0]++;
            nPartners[idx1][1]++;
        }
    }

    public int basisSize(){
    	return 1;
    }

    public boolean haveTarget(IAtom a){
        if (parentGroup.getChildList().getAtom(a.getIndex()) != a) return false;
        if (a.getIndex() >= partners.length) return false;
        int[][] aPartners = partners[a.getIndex()];
        if (aPartners[0].length + aPartners[1].length == 0) return false;
    	return true;
    }

    public void setTarget(IAtom a){
    	target = a;
    	unset();
    }

    public void setBasis(IMoleculeList parent) {
        if (parent == null) {
            parentGroup = null;
        }
        else {
            parentGroup = parent.getMolecule(0);
        }
        unset();
    }

    public int size() {
        if (target == null) return index.length;
        if (!haveTarget(target)) return 0;
        return partners[target.getIndex()][0].length + partners[target.getIndex()][1].length;
    }

    /**
     * Sets iterator to a state where hasNext() returns false.
     */
    public void unset() {
        cursor = index.length;
    }

    /**
     * Resets iterator to a state where hasNext is true, if atoms in pair are
     * not null.
     */
    public void reset() {
        cursor = index.length;
        if (parentGroup != null) {
            cursor = 0;
            if (target!=null && !haveTarget(target)) cursor = index.length;
        }
        
    }

    /**
     * Returns the iterator's pair and unsets iterator.
     */
    public IAtomList next() {
        if (target != null){
            int[][] tPartners = partners[target.getIndex()];
            int t0Length = tPartners[0].length;
            if (cursor < t0Length) {
                pair.atom0 = target;
                pair.atom1 = parentGroup.getChildList().getAtom(tPartners[0][cursor]);
                cursor++;
                return pair;
            }
            if (cursor < t0Length + tPartners[1].length) {
                pair.atom1 = target;
                pair.atom0 = parentGroup.getChildList().getAtom(tPartners[1][cursor-t0Length]);
                cursor++;
                return pair;
            }
            return null;
        }

        if (cursor >= index.length){
            return null;
        }
        pair.atom0 = parentGroup.getChildList().getAtom(index[cursor][0]);
        pair.atom1 = parentGroup.getChildList().getAtom(index[cursor][1]);
        cursor++;
        return pair;
    }

    /**
     * Returns 2, indicating that this is a pair iterator.
     */
    public final int nBody() {
        return 2;
    }

    protected final int [][] index;
    protected final int [][][] partners;
    protected int cursor;
    protected IMolecule parentGroup;
    protected IAtom target;
    protected final AtomPair pair;
}

