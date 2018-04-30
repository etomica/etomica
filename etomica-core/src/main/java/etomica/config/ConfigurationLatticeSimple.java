/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.config;

import etomica.action.MoleculeActionTranslateTo;
import etomica.box.Box;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.IndexIteratorRectangular;
import etomica.lattice.IndexIteratorSizable;
import etomica.lattice.SpaceLattice;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculePositionGeometricCenter;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Constructs configuration that has the molecules placed on the sites of a
 * lattice. Molecules are placed sequentially on the sites of the lattice in an
 * order specified by an index iterator set at construction.
 * <p>
 * Iteration over the indexes yields integer arrays, and with each iteration the
 * array is passed to the site method of the lattice which returns the position
 * for placement of the next molecule. Each array index is iterated to a maximum
 * value determined by the number of molecules to be placed, the dimensions of
 * the box in which they are placed, and the lattice constants of the lattice.
 * <p>
 * An instance of this class may be configured to place atoms such that they
 * uniformly fill the volume of the box. It will attempt this by scaling the
 * lattice constants of the configuration in an appropriate way. Success in
 * getting a good spatial distribution may vary.
 * <p>
 * An instance can also be configured to remember the indices used to get
 * lattice position for each molecule that is placed. This can be useful if it
 * is desired to associate each molecule with a lattice site.
 */
public class ConfigurationLatticeSimple implements Configuration, java.io.Serializable {

    /**
     * Constructs class using instance of IndexIteratorRectangular as the default
     * index iterator.
     */
    public ConfigurationLatticeSimple(SpaceLattice lattice, Space space) {
        this(lattice, new IndexIteratorRectangular(lattice.D()), space);
    }

    /**
     * Construct class that will place atoms on sites of the given lattice,
     * proceeding in the order resulting from iteration through the given index
     * iterator.
     */
    public ConfigurationLatticeSimple(SpaceLattice lattice,
            IndexIteratorSizable indexIterator, Space space) {
        if(indexIterator.getD() != lattice.D()) {
            throw new IllegalArgumentException("Dimension of index iterator and lattice are incompatible");
        }
        this.space = space;
        this.lattice = lattice;
        this.indexIterator = indexIterator;
        atomActionTranslateTo = new MoleculeActionTranslateTo(lattice.getSpace());
        atomActionTranslateTo.setAtomPositionDefinition(new MoleculePositionGeometricCenter(space));
    }

    /**
     * Places the molecules in the given box on the positions of the
     * lattice.  
     */
    public void initializeCoordinates(Box box) {
        IMoleculeList moleculeList = box.getMoleculeList();
        int sumOfMolecules = moleculeList.size();
        if (sumOfMolecules == 0) {
            return;
        }
        int basisSize = 1;
        if (lattice instanceof BravaisLatticeCrystal) {
            basisSize = ((BravaisLatticeCrystal)lattice).getBasis().getScaledCoordinates().length;
        }
        int nCells = (int) Math.ceil((double) sumOfMolecules
                / (double) basisSize);

        // determine scaled shape of simulation volume
        Vector dim = space.makeVector();
        Boundary boundary = box.getBoundary();
        for (int i=0; i<space.D(); i++) {
            Vector edgeVector = boundary.getEdgeVector(i);
            dim.setX(i,Math.sqrt(edgeVector.squared()));
        }
        Vector shape = space.makeVector();
        shape.E(dim);
        Vector latticeConstantV = Vector.of(lattice.getLatticeConstants());
        shape.DE(latticeConstantV);

        // determine number of cells in each direction
        int[] latticeDimensions = calculateLatticeDimensions(nCells, shape);
        if (indexIterator.getD() > latticeDimensions.length) {
            int[] iteratorDimensions = new int[latticeDimensions.length+1];
            System.arraycopy(latticeDimensions, 0, iteratorDimensions, 0,
                    latticeDimensions.length);
            iteratorDimensions[latticeDimensions.length] = basisSize;
            indexIterator.setSize(iteratorDimensions);
        }
        else {
            indexIterator.setSize(latticeDimensions);
        }

        // Place molecules
        indexIterator.reset();
        Vector offset = space.makeVector();
        offset.Ea1Tv1(-0.5, box.getBoundary().getBoxSize());
        Vector destinationVector = atomActionTranslateTo.getDestination();
        int nMolecules = moleculeList.size();
        for (int i=0; i<nMolecules; i++) {
            IMolecule a = moleculeList.get(i);
            // initialize coordinates of child atoms
            a.getType().initializeConformation(a);

            destinationVector.Ev1Pv2((Vector)lattice.site(indexIterator.next()), offset);
            atomActionTranslateTo.actionPerformed(a);
        }
    }

    protected int[] calculateLatticeDimensions(int nCells, Vector shape) {
        int dimLeft = shape.getD();
        int nCellsLeft = nCells;
        int[] latticeDimensions = new int[shape.getD()];
        while (dimLeft > 0) {
            double smin = Double.POSITIVE_INFINITY;
            int dmin = 0;
            double product = 1.0;
            for (int idim = 0; idim < shape.getD(); idim++) {
                if (latticeDimensions[idim] > 0)
                    continue;
                if (shape.getX(idim) < smin) {
                    smin = shape.getX(idim);
                    dmin = idim;
                }
                product *= shape.getX(idim);
            }
            // round off except for last dimension (then round up)
            if (dimLeft > 1) {
                latticeDimensions[dmin] = (int) Math.round(shape.getX(dmin)
                        * Math.pow((nCellsLeft / product), 1.0 / dimLeft));
                if (latticeDimensions[dmin] == 0){
                	latticeDimensions[dmin] = 1;
                }
            } else {
                latticeDimensions[dmin] = (int) Math.ceil(shape.getX(dmin)
                        * nCellsLeft / product);
            }
            nCellsLeft = (nCellsLeft + latticeDimensions[dmin] - 1)
                    / latticeDimensions[dmin];
            dimLeft--;
        }
        return latticeDimensions;
    }

    protected final SpaceLattice lattice;
    protected final IndexIteratorSizable indexIterator;
    protected final MoleculeActionTranslateTo atomActionTranslateTo;
    private final Space space;
    private static final long serialVersionUID = 2L;
}
