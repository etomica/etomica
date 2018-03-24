/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.lattice;

import etomica.potential.IteratorDirective;

/**
 * An abstract iterator that generates the neighboring sites of a given site.
 * Subclasses provide the specification of the "neighbors" of a site by a suitable
 * implementation of the updateNeighborList method.
 * <p>
 * Iterator can be configured to yield the up- and/or down-neighbors of the site.
 * The definition of "up" and "down" neighbors is made by the subclass, within
 * the restriction that if site B is an up-neighbor of site A, then A will
 * be a down-neighbor of site B.  Note that because of the way periodic boundaries
 * are treated, a transitive relation is not guaranteed,
 * meaning that (B upNbr of A) and (C upNbr of B) does not imply (C upNbr of A).
 * <p>
 * The central site is not included among the iterates (not a self-neighbor).
 */
public abstract class RectangularLatticeNbrIterator implements SiteIterator {

    protected final int D;
    protected final boolean[] isPeriodic;
    final int[] centralSite;
    private final int[] latticeIndex;
    protected boolean needNeighborUpdate = true;
    protected int[] neighbors;
    protected RectangularLattice lattice;
    protected int cursor;
    protected int neighborCount;
    protected boolean doUp, doDown;
    private IteratorDirective.Direction direction;

    /**
     * Constructs iterator that needs to be configured before use.  Must specify
     * the central site and the range before reset and iteration.
     */
    public RectangularLatticeNbrIterator(int D) {
        super();
        this.D = D;
        centralSite = new int[D];
        cursor = Integer.MAX_VALUE;
        isPeriodic = new boolean[D];
        for (int i = 0; i < D; i++) {
            isPeriodic[i] = true;
        }
        latticeIndex = new int[D];
    }

    protected abstract void updateNeighborList();

    public final int D() {
        return D;
    }

    /**
     * Identifies the site whose neighbors will be given on iteration.
     * Requires subsequent call to reset before iteration.
     *
     * @param index the coordinate of the central site
     */
    public void setSite(int[] index) {
        if (index.length != D) throw new IllegalArgumentException("Incorrect length of array passed to setSite");
        for (int i = D - 1; i >= 0; i--) centralSite[i] = index[i];
        needNeighborUpdate = true;
        unset();
    }

    /**
     * Indicates whether the iterator has another site.
     */
    public final boolean hasNext() {
        return cursor < neighborCount;
    }

    /**
     * Resets the iterator to loop through its iterates again.  This
     * must be done after any call to setRange, setDirection, or setSite.
     */
    //other iterators assume that this call is inexpensive if neighbor list is not being updated
    public void reset() {
        if (lattice == null) {
            cursor = 0;
            neighborCount = 0;
            return;
        }
        if (needNeighborUpdate) updateNeighborList();
        cursor = 0;
    }

    /**
     * Puts iterator in a state in which hasNext() returns false.
     */
    public void unset() {
        cursor = neighborCount;
    }

    /**
     * Returns the index of the next iterate while advancing the
     * iterator.
     */
    public int[] nextIndex() {
        if (!hasNext()) {
            return null;
        }
        lattice.latticeIndex(neighbors[cursor++], latticeIndex);
        return latticeIndex;
    }

    /**
     * Returns the next site in the iteration sequence, or null
     * if hasNext is false.
     */
    public Object next() {
//            if(!hasNext()) return null;
        //           currentPbc = pbc[cursor];
//            return lattice.sites[neighbors[cursor++]];
        return hasNext() ? lattice.sites()[neighbors[cursor++]] : null;
    }

    /**
     * Returns the next iterate without advancing the iterator.
     */
    public Object peek() {
        return hasNext() ? lattice.sites()[neighbors[cursor]] : null;
    }

    /**
     * The number of iterates returned by this iterator in its current state.
     */
    public int size() {
        return neighborCount;
    }

    /**
     * Gets the direction (up or down list) of the neighbors given by the iterator.
     * A null value indicates that all neighbors will be returned.
     */
    public IteratorDirective.Direction getDirection() {
        return direction;
    }

    /**
     * Sets the iterator to return only up-list neighbors or down-list neighbors.
     * A null value indicates that all neighbors are to be returned.
     */
    public void setDirection(IteratorDirective.Direction direction) {
        this.direction = direction;
        doUp = (direction != IteratorDirective.Direction.DOWN);//also handles case where
        doDown = (direction != IteratorDirective.Direction.UP);//direction is null
        needNeighborUpdate = true;
        unset();
    }

    /**
     * Returns the lattice from which the iterates are given.
     */
    public RectangularLattice getLattice() {
        return lattice;
    }

    /**
     * Sets the lattice from which the iterates are given.  Dimension
     * of lattice must be the same as that specified to constructor. Also,
     * size of lattice must be compatible with range of neighbor interactions.
     */
    public void setLattice(FiniteLattice lattice) {
        this.lattice = (RectangularLattice) lattice;
        if (lattice != null && lattice.D() != D)
            throw new IllegalArgumentException("Iterator given lattice with incompatible dimension");
        needNeighborUpdate = true;
    }

    public void setPeriodicity(boolean[] periodicity) {
        for (int i = 0; i < D; i++) {
            isPeriodic[i] = periodicity[i];
        }
    }
}
