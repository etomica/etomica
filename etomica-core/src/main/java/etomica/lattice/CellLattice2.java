package etomica.lattice;

import etomica.nbr.cell.Cell;
import etomica.potential.IteratorDirective;
import etomica.space.Space;
import etomica.space.Vector;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class CellLattice2 implements FiniteLattice<Cell> {

    private final int d;
    private final int[] jumpCount; // jumpCount[i] gives the number of sites skipped when the i-th index is incremented by 1
    private final int[] size; // number of cells in each dimension
    private final double[] cellSize; // size of a cell in each dimension
    private final Vector dimensions;
    private Cell[] sites;
    private int[][] upNeighbors; // for each cell index, list the indices of its neighbors

    public CellLattice2(Space space, Vector dimensions, int[] size, double neighborRange) {
        this.d = space.getD();
        this.jumpCount = new int[d];
        jumpCount[d - 1] = 1;
        this.size = new int[d];
        this.cellSize = new double[d];
        this.dimensions = space.makeVector();
        this.dimensions.E(dimensions);
        this.setSize(size); // initializes sites
        this.upNeighbors = new int[sites.length][];

        // compute up neighbors
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

    public int[][] getSiteUpNeighbors() {
        return upNeighbors;
    }

    @Override
    public Cell[] sites() {
        return this.sites;
    }

    @Override
    public int[] getSize() {
        return this.size;
    }

    @Override
    public void setSize(int[] size) {
        System.arraycopy(size, 0, this.size, 0, d);
        for (int i = d - 1; i > 0; i--) {
            jumpCount[i - 1] = jumpCount[i] * this.size[i];
        }

        this.sites = new Cell[jumpCount[0] * this.size[0]];
        for (int i = 0; i < sites.length; i++) {
            sites[i] = new Cell(i);
        }

    }

    @Override
    public int D() {
        return this.d;
    }

    /**
     * Returns the cell in which the given point lies. An ArrayIndexOutOfBounds
     * exception is thrown if the point lies outside the bounds of the lattice
     * (less than zero or greater than dimensions). The dimension of the vector
     * is assumed to be consistent with the dimension D of the lattice (i.e.,
     * r.D() == this.D()) but this is not checked.
     */
    public Object site(Vector r) {
        int idx1D = 0;
        for (int i = 0; i < d; i++) {
            int j = ((int) (size[i] * (r.getX(i) / dimensions.getX(i) + 0.5)));
            if (j == -1) j = 0;
            else if (j == size[i]) j = size[i] - 1;
            idx1D += j * jumpCount[i];
        }
        return sites[idx1D];
    }

    @Override
    public Cell site(int[] index) {
        return sites[arrayIndex(index)];
    }

    public final int arrayIndex(int[] latticeCoords) {
        int idx = 0;
        for (int i = 0; i < d; i++) {
            idx += latticeCoords[i] * jumpCount[i];
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
}
