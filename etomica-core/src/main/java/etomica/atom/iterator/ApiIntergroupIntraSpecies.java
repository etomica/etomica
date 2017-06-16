/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;

/**
 * Iterator that returns pairs formed using two different basis atoms, so that
 * the iterates are taken from two different groups.  This subclass of
 * ApiIntergroup is specifically designed to return atoms of different types
 * from different molecules of the same species.  For methane, it
 * would return both C-H and H-C.
 */
public class ApiIntergroupIntraSpecies extends ApiIntergroup {

    /**
     * Constructs a pair iterator that returns iterates from the given
     * pairIterator, which is expected to contain two basis-dependent 
     * iterators.
     */
    public ApiIntergroupIntraSpecies(AtomIteratorBasisDependent outer, AtomIteratorBasisDependent inner) {
        super(outer, inner);
    }

    /**
     * Returns the next pair of atoms. The same AtomPair instance is returned
     * every time, but the Atoms it holds are (of course) different for each
     * iterate.
     */
    public IAtomList next() {
        //Advance the inner loop, if it is not at its end.
        IAtom nextInner = aiInner.nextAtom();
        while (nextInner != null && nextInner.getType() == pair.atom0.getType()) {
            nextInner = aiInner.nextAtom();
        }
        if (nextInner != null) {
            pair.atom1 = nextInner;
            return pair;
        }

        while (true) {
            //Advance the outer loop, if the inner loop has reached its end.
            IAtom nextOuter = aiOuter.nextAtom();
            if (nextOuter == null) {
                return null;
            }
            pair.atom0 = nextOuter;

            aiInner.reset();
            while (nextInner != null && nextInner.getType() == nextOuter.getType()) {
                nextInner = aiInner.nextAtom();
            }

            if (nextInner != null) {
                pair.atom1 = nextInner;
                return pair;
            }
        }
    }
}
