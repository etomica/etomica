/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.lattice;


import etomica.potential.IteratorDirective;

import java.util.ArrayList;
import java.util.List;

/**
 * Basic implementation of the AbstractLattice interface, providing construction
 * and access of sites for a lattice of arbitrary dimension. Lattice is
 * retangular in the sense that the size in one dimension does not depend on the
 * index in another dimension (e.g., it cannot be triangular). <p>
 * RectangularLatticeNbrIterator defines a configurable neighbor iterator that
 * returns the sites it defines as neighbors of a given site.
 *
 * @see RectangularLatticeNbrIterator
 */

//Internally, sites are stored in a 1-D array of objects, and are accessed by
//unrolling the index specification to determine the storage-array index.

// Example showing internal ordering of elements
//  0     1     2     3     4     5     6     7     8     9    10    11   arrayIndex
//(000) (001) (002) (010) (011) (012) (100) (101) (102) (110) (111) (112) latticeIndex
//  for this example, size = {2, 2, 3}, jumpCount = {6, 3, 1}
//  note that number of sites = size[0]*jumpCount[0]

public class RectangularLattice implements FiniteLattice {

    protected final int[] size;
    //  jumpCount[i] gives the number of sites skipped when the i-th index is incremented by 1
    protected final int[] jumpCount;
    protected final int d;
    protected Object[] sites;
    protected SiteFactory siteFactory;

    private double neighborRange = 1.0;
    private int[][] upNeighbors;

    /**
     * Constructs a lattice of the given dimension (D) with sites
     * made from the given factory.  Upon construction, lattice contains
     * no sites; these are created when setSize is invoked.
     */
    public RectangularLattice(int D, SiteFactory siteFactory) {
        this.d = D;
        jumpCount = new int[D];
        jumpCount[D - 1] = 1;
        size = new int[D];
        this.siteFactory = siteFactory;
        //do not create lattice with default size because siteFactory  might not yet be ready
    }

    private void computeUpNeighbors() {
        this.upNeighbors = new int[sites.length][];
        CellLattice.NeighborIterator iter = new CellLattice.NeighborIterator(d, neighborRange);
        iter.setLattice(this);
        iter.setDirection(IteratorDirective.Direction.UP);
        for (int i = 0; i < sites().length; i++) {
            iter.setSite(this.latticeIndex(i));
            iter.reset();
            List<Integer> nbrIndices = new ArrayList<>();
            while(iter.hasNext()) {
                nbrIndices.add(this.arrayIndex(iter.nextIndex()));
            }

            upNeighbors[i] = new int[nbrIndices.size()];
            for (int j = 0; j < nbrIndices.size(); j++) {
                upNeighbors[i][j] = nbrIndices.get(j);
            }
        }
    }

    public int[][] getUpNeighbors() {
        return upNeighbors;
    }

    @Override
    public final int D() {
        return d;
    }

    @Override
    public Object[] sites() {
        return sites;
    }

    /**
     * Returns the instance of the object associated with the given index.
     * Repeated calls with the same index will return the same instance, and
     * calls with different indexes will return different instances.
     */
    @Override
    public Object site(int[] index) {
        return sites[arrayIndex(index)];
    }

    /**
     * Returns the index in the 1-d array for the site corresponding
     * to the given lattice index.
     */
    public final int arrayIndex(int[] index) {
        int idx = 0;
        for (int i = 0; i < d; i++) {
            idx += index[i] * jumpCount[i];
        }
        return idx;
    }

    /**
     * Returns the lattice index given the 1-d array index; reverses
     * the effect of arrayIndex method.
     */
    public int[] latticeIndex(int index) {
        int[] latticeIndex = new int[d];
        latticeIndex(index, latticeIndex);
        return latticeIndex;
    }

    public void latticeIndex(int index, int[] latticeIndex) {
        for (int i = 0; i < d; i++) {
            latticeIndex[i] = index / jumpCount[i];
            index -= latticeIndex[i] * jumpCount[i];
        }
    }

    /* (non-Javadoc)
     * @see etomica.lattice.AbstractLattice#getDimensions()
     */
    @Override
    public final int[] getSize() {
        return size;
    }

    /**
     * Sets the number lattice size (specified via the largest index in each dimension), and
     * rebuilds the all sites using the site factory.
     *
     * @param newSize array giving the number of index values in each dimension
     * @throws IllegalArgumentException if length of given array is not equal to D
     */
    @Override
    public void setSize(int[] newSize) {
        if (newSize.length != d) throw new IllegalArgumentException("Incorrect dimension dimension");
        System.arraycopy(newSize, 0, size, 0, d);
        for (int i = d - 1; i > 0; i--) {
            jumpCount[i - 1] = jumpCount[i] * size[i];
        }
        sites = new Object[jumpCount[0] * size[0]];
        int[] idx = new int[d];
        idx[d - 1] = -1;
        for (int i = 0; i < sites.length; i++) {
            increment(idx);
            sites[i] = siteFactory.makeSite(this, idx);
        }

        this.computeUpNeighbors();
    }

    //method used by setDimensions method to cycle the index array through its values
    protected void increment(int[] idx) {
        int d = this.d - 1;
        idx[d]++;
        while (idx[d] == size[d] && d > 0) {//replaces recursive call
            idx[d] = 0;
            idx[--d]++;//decrement d, then increment idx
        }
    }

    public void setNeighborRange(double neighborRange) {
        this.neighborRange = neighborRange;
        this.computeUpNeighbors();
    }

    /**
     * Iterates over all sites of the lattice.
     */
    public static class Iterator implements SiteIterator, java.io.Serializable {

        private static final long serialVersionUID = 1L;
        private final int[] idx;//index of the most recently returned iterate
        private int cursor = Integer.MAX_VALUE;
        private RectangularLattice lattice;
        private int size = 0;

        public Iterator(int D) {
            idx = new int[D];
        }

        @Override
        public boolean hasNext() {
            return cursor < size;
        }

        @Override
        public Object next() {
            if (hasNext()) {
                lattice.increment(idx);
                return lattice.sites[cursor++];
            }
            return null;
        }

        @Override
        public int[] nextIndex() {
            if (hasNext()) {
                lattice.increment(idx);
                cursor++;
                return idx;
            }
            return null;
        }

        @Override
        public Object peek() {
            return hasNext() ? lattice.sites[cursor] : null;
        }

        @Override
        public void reset() {
            size = size();
            cursor = 0;
            for (int i = 0; i < idx.length; i++) idx[i] = 0;
            idx[idx.length - 1] = -1;
        }

        @Override
        public int size() {
            return (lattice != null) ? lattice.sites.length : 0;
        }

        @Override
        public void unset() {
            cursor = Integer.MAX_VALUE;
        }

        @Override
        public void setLattice(FiniteLattice lattice) {
            this.lattice = (RectangularLattice) lattice;
            unset();
        }
    }
}
