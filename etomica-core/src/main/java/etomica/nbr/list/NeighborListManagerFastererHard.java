/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.nbr.list;

import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.box.Box;
import etomica.potential.BondingInfo;
import etomica.potential.IPotentialAtomic;
import etomica.potential.IPotentialHard;
import etomica.potential.compute.NeighborIterator;
import etomica.potential.compute.NeighborManagerHard;
import etomica.space.Vector;
import etomica.species.SpeciesManager;

/**
 * Neighbor manager that that uses neighbor lists and also stores the state of
 * neighbor pairs for hard MD.
 */
public class NeighborListManagerFastererHard extends NeighborListManagerFasterer implements NeighborManagerHard {

    public int[][] nbrState;
    private final NeighborIteratorListHard neighborIterator;

    public NeighborListManagerFastererHard(SpeciesManager sm, Box box, int cellRange, double nbrRange, BondingInfo bondingInfo) {
        super(sm, box, cellRange, nbrRange, bondingInfo);
        setDoDownNeighbors(true);
        neighborIterator = new NeighborIteratorListHard(this, box);
    }

    @Override
    public NeighborIterator makeNeighborIterator() {
        return this.neighborIterator;
    }

    @Override
    protected void realloc() {
        super.realloc();
        nbrState = new int[box.getLeafList().size()][maxNab];
    }

    @Override
    protected void newDownNeighbor(int j, int i, int upSlot, int downSlot) {
        super.newDownNeighbor(j, i, upSlot, downSlot);
        nbrState[j][downSlot] = nbrState[i][upSlot];
    }

    @Override
    protected int checkNbrPair(int i, int j, IAtom iAtom, IAtom jAtom, double rc2, Vector jbo, IPotentialAtomic[] iPotentials) {
        // we need to override this whole method so we can grab and store the pair state for neighbors
        if (iPotentials[jAtom.getType().getIndex()] == null) return 0;

        if (bondingInfo.skipBondedPair(isPureAtoms, iAtom, jAtom)) return 0;

        Vector dr = space.makeVector();
        Vector ri = iAtom.getPosition();
        Vector rj = jAtom.getPosition();
        dr.Ev1Mv2(rj, ri);
        dr.PE(jbo);
        double r2 = dr.squared();
        if (r2 > rc2) return 0;
        if (numAtomNbrsUp[i] >= maxNab) return 1;
        nbrs[i][numAtomNbrsUp[i]] = j;
        int state = ((IPotentialHard) iPotentials[jAtom.getType().getIndex()]).getState((IAtomKinetic) iAtom, (IAtomKinetic) jAtom, dr);
        nbrState[i][numAtomNbrsUp[i]] = state;
        nbrBoxOffsets[i][numAtomNbrsUp[i]] = jbo;
        numAtomNbrsUp[i]++;
        return 0;
    }

    @Override
    public void setPairState(int i, int j, int state) {
        // a collision happened
        for (int ii = 0; ii < numAtomNbrsUp[i]; ii++) {
            if (nbrs[i][ii] == j) {
                nbrState[i][ii] = state;
                break;
            }
        }
        for (int jj = maxNab - 1; jj > maxNab - 1 - numAtomNbrsDn[j]; jj--) {
            if (nbrs[j][jj] == i) {
                nbrState[j][jj] = state;
                break;
            }
        }
    }
}
