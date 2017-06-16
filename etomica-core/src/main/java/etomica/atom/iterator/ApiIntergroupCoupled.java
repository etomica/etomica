/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.atom.IAtomList;

/**
 * Iterator that returns coupled iterates; the first pair contains the first
 * atom from each basis.  The second pair has the second atom from each basis.
 *
 * @author Andrew Schultz
 */
public class ApiIntergroupCoupled extends ApiIntergroup {

    public ApiIntergroupCoupled() {
    }

    public ApiIntergroupCoupled(AtomIteratorBasis outer, AtomIteratorBasis inner) {
        super(outer, inner);
    }

    /**
     * Returns the number of pairs given by this iterator. Not dependent on
     * state of hasNext.
     */
    public int size() {
        return aiOuter.size();
    }

    /**
     * Returns the next pair of atoms. The same AtomPair instance is returned
     * every time, but the Atoms it holds are (of course) different for each
     * iterate.
     */
    public IAtomList next() {
        //Advance the inner loop, if it is not at its end.
        pair.atom0 = aiOuter.nextAtom();
        if (pair.atom0 == null) {
            return null;
        }
        pair.atom1 = aiInner.nextAtom();
        return pair;
    }

    public void reset() {
        aiOuter.reset();
        aiInner.reset();
    }
}
