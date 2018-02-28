/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.lattice;


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

public class RectangularLattice implements FiniteLattice, java.io.Serializable {

    private static final long serialVersionUID = 1L;
    protected final int[] size;
    //  jumpCount[i] gives the number of sites skipped when the i-th index is incremented by 1
    protected final int[] jumpCount;
    protected final int D;
    protected Object[] sites;
    protected SiteFactory siteFactory;

    /**
     * Constructs a lattice of the given dimension (D) with sites
     * made from the given factory.  Upon construction, lattice contains
     * no sites; these are created when setSize is invoked.
     */
    public RectangularLattice(int D, SiteFactory siteFactory) {
        this.D = D;
        jumpCount = new int[D];
        jumpCount[D - 1] = 1;
        size = new int[D];
        this.siteFactory = siteFactory;
        //do not create lattice with default size because siteFactory  might not yet be ready
    }

    /* (non-Javadoc)
     * @see etomica.lattice.AbstractLattice#D()
     */
    public final int D() {
        return D;
    }

    /* (non-Javadoc)
     * @see etomica.lattice.AbstractLattice#siteList()
     */
    public Object[] sites() {
        return sites;
    }

    /**
     * Returns the instance of the object associated with the given index.
     * Repeated calls with the same index will return the same instance, and
     * calls with different indexes will return different instances.
     */
    public Object site(int[] index) {
        return sites[arrayIndex(index)];
    }

    /**
     * Returns the index in the 1-d array for the site corresponding
     * to the given lattice index.
     */
    public final int arrayIndex(int[] index) {
        int idx = 0;
        for (int i = 0; i < D; i++) {
            idx += index[i] * jumpCount[i];
        }
        return idx;
    }

    /**
     * Returns the lattice index given the 1-d array index; reverses
     * the effect of arrayIndex method.
     */
    public int[] latticeIndex(int index) {
        int[] latticeIndex = new int[D];
        latticeIndex(index, latticeIndex);
        return latticeIndex;
    }

    public void latticeIndex(int index, int[] latticeIndex) {
        for (int i = 0; i < D; i++) {
            latticeIndex[i] = index / jumpCount[i];
            index -= latticeIndex[i] * jumpCount[i];
        }
    }

    /* (non-Javadoc)
     * @see etomica.lattice.AbstractLattice#getDimensions()
     */
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
    public void setSize(int[] newSize) {
        if (newSize.length != D) throw new IllegalArgumentException("Incorrect dimension dimension");
        System.arraycopy(newSize, 0, size, 0, D);
        for (int i = D - 1; i > 0; i--) {
            jumpCount[i - 1] = jumpCount[i] * size[i];
        }
        sites = new Object[jumpCount[0] * size[0]];
        int[] idx = new int[D];
        idx[D - 1] = -1;
        for (int i = 0; i < sites.length; i++) {
            increment(idx);
            sites[i] = siteFactory.makeSite(this, idx);
        }
    }

    //method used by setDimensions method to cycle the index array through its values
    protected void increment(int[] idx) {
        int d = D - 1;
        idx[d]++;
        while (idx[d] == size[d] && d > 0) {//replaces recursive call
            idx[d] = 0;
            idx[--d]++;//decrement d, then increment idx
        }
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

        public boolean hasNext() {
            return cursor < size;
        }

        public Object next() {
            if (hasNext()) {
                lattice.increment(idx);
                return lattice.sites[cursor++];
            }
            return null;
        }

        public int[] nextIndex() {
            if (hasNext()) {
                lattice.increment(idx);
                cursor++;
                return idx;
            }
            return null;
        }

        public Object peek() {
            return hasNext() ? lattice.sites[cursor] : null;
        }

        public void reset() {
            size = size();
            cursor = 0;
            for (int i = 0; i < idx.length; i++) idx[i] = 0;
            idx[idx.length - 1] = -1;
        }

        public int size() {
            return (lattice != null) ? lattice.sites.length : 0;
        }

        public void unset() {
            cursor = Integer.MAX_VALUE;
        }

        public void setLattice(FiniteLattice lattice) {
            this.lattice = (RectangularLattice) lattice;
            unset();
        }
    }
}
