/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import java.util.ArrayList;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.api.IMoleculeList;
import etomica.atom.AtomArrayList;

/**
 * Iterator that returns coupled iterates of any size; the first set contains
 * the first atom from each basis.  The second set has the second atom from
 * each basis.
 * 
 * @author Andrew Schultz
 */
public class ANIntergroupCoupled implements AtomsetIteratorBasisDependent {

    public ANIntergroupCoupled(int nBody) {
        atoms = new AtomArrayList(nBody);
        atomLists = new ArrayList<IAtomList>(nBody);
        this.nBody = nBody;
        unset();
    }

    public void unset() {
        counter = -1;
    }

    public int size() {
        return nLeafAtoms;
    }

    public final int nBody() {
        return nBody;
    }

    public IAtomList next() {
        if (counter >= nLeafAtoms) {
            return null;
        }
        atoms.clear();
        for (int i=0; i<nBody; i++) {
            atoms.add(atomLists.get(i).getAtom(counter));
        }
        counter++;
        return atoms;
    }

    public void setTarget(IAtom newTargetAtom) {
        // we are basis dependent because we have setBasis, but this
        // method should not be called.
        throw new RuntimeException("nope");
    }

    public boolean haveTarget(IAtom target) {
        // we are basis dependent because we have setBasis, but this
        // method should not be called.
        throw new RuntimeException("nope");
    }

    /**
     * Specifies the basis, which identifies the molecules subject to
     * iteration.
     */
    public void setBasis(IMoleculeList basisMolecules) {
        if (basisMolecules.getMoleculeCount() != nBody) {
            throw new IllegalArgumentException();
        }
        atomLists.clear();
        for (int i=0; i<nBody; i++) {
            atomLists.add(basisMolecules.getMolecule(i).getChildList());
        }
        nLeafAtoms = atomLists.get(0).getAtomCount();
    }

    public void reset() {
        counter = 0;
    }

    /**
     * Returns 2, indicating that the setBasis method expects an array of two
     * atoms.
     */
    public int basisSize() {
        return nBody;
    }

    private static final long serialVersionUID = 1L;
    protected final AtomArrayList atoms;
    protected int counter, nLeafAtoms;
    protected final ArrayList<IAtomList> atomLists;
    protected boolean needSetupIterators = true;
    protected final int nBody;
}
