/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.paracetamol;

import etomica.action.MoleculeActionTranslateTo;
import etomica.action.MoleculeChildAtomAction;
import etomica.space.Vector;
import etomica.box.Box;
import etomica.atom.IMolecule;
import etomica.atom.IMoleculeList;
import etomica.config.Configuration;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.IndexIteratorRectangular;
import etomica.lattice.IndexIteratorSizable;
import etomica.lattice.SpaceLattice;
import etomica.space.Space;
import etomica.space.Tensor;

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
public class ConfigurationMonoclinicLattice implements Configuration, java.io.Serializable {

	private final static String APP_NAME = "Configuration Monoclinic Lattice";
	private final Space space;

    /**
     * Constructs class using instance of IndexIteratorRectangular as the default
     * index iterator.
     */
    public ConfigurationMonoclinicLattice(SpaceLattice lattice, Space _space) {
        this(lattice, new IndexIteratorRectangular(lattice.D()), _space);
    }

    /**
     * Construct class that will place atoms on sites of the given lattice,
     * proceeding in the order resulting from iteration through the given index
     * iterator.
     */
    public ConfigurationMonoclinicLattice(SpaceLattice lattice,
            IndexIteratorSizable indexIterator, Space _space) {
        if(indexIterator.getD() != lattice.D()) {
            throw new IllegalArgumentException("Dimension of index iterator and lattice are incompatible");
        }
        this.lattice = lattice;
        this.indexIterator = indexIterator;
        this.space = _space;
        atomActionTranslateTo = new MoleculeActionTranslateTo(lattice.getSpace());
        atomGroupAction = new MoleculeChildAtomAction(new AtomActionTransformed(lattice.getSpace()));
    }

    /**
     * Places the molecules in the given box on the positions of the
     * lattice.  
     */
    public void initializeCoordinates(Box box) {
        IMoleculeList moleculeList = box.getMoleculeList();
        int sumOfMolecules = moleculeList.getMoleculeCount();
        if (sumOfMolecules == 0) {
            return;
        }
        int basisSize = 4;
        
        if (lattice instanceof BravaisLatticeCrystal) {
            basisSize = ((BravaisLatticeCrystal)lattice).getBasis().getScaledCoordinates().length;
        }
        int nCells = (int) Math.ceil((double) sumOfMolecules
                / (double) basisSize);

        // determine scaled shape of simulation volume
        Vector shape = space.makeVector();
        shape.E(box.getBoundary().getBoxSize());
        Vector latticeConstantV = space.makeVector(lattice.getLatticeConstants());
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

        // determine lattice constant
        Vector latticeScaling = space.makeVector();
        if (rescalingToFitVolume) {
            // in favorable situations, this should be approximately equal
            // to 1.0
            latticeScaling.E(box.getBoundary().getBoxSize());
            latticeScaling.DE(latticeConstantV);
            latticeScaling.DE(space.makeVector(latticeDimensions));
        } else {
            latticeScaling.E(1.0);
        }

        // determine amount to shift lattice so it is centered in volume
        Vector offset = space.makeVector();
        offset.E(box.getBoundary().getBoxSize());
        Vector vectorOfMax = space.makeVector();
        Vector vectorOfMin = space.makeVector();
        Vector site = space.makeVector();
        vectorOfMax.E(Double.NEGATIVE_INFINITY);
        vectorOfMin.E(Double.POSITIVE_INFINITY);

        // XXX this can do strange things. it's probably not needed for 
        // periodic boundaries, but gets the atoms off the boundaries for 
        // non-periodic boundaries
        indexIterator.reset();

        while (indexIterator.hasNext()) {
            site.E((Vector) lattice.site(indexIterator.next()));
            site.TE(latticeScaling);
            for (int i=0; i<site.getD(); i++) {
                vectorOfMax.setX(i, Math.max(site.getX(i),vectorOfMax.getX(i)));
                vectorOfMin.setX(i, Math.min(site.getX(i),vectorOfMin.getX(i)));
            }
        }
        offset.Ev1Mv2(vectorOfMax, vectorOfMin);
        offset.TE(-0.5);
        offset.ME(vectorOfMin);

        myLat = new MyLattice(lattice, latticeScaling, offset);

        // Place molecules
        indexIterator.reset();

    	ConformationParacetamolMonoclinic regConfig = new ConformationParacetamolMonoclinic(lattice.getSpace());
    	Vector cellPosition = null;
    	Tensor t = lattice.getSpace().makeTensor();

    	for (int iMolecule = 0; iMolecule<moleculeList.getMoleculeCount(); iMolecule++) {
    	    IMolecule molecule = moleculeList.getMolecule(iMolecule);
    		    
            int[] ii = indexIterator.next();
            
            regConfig.initializePositions(molecule.getChildList());
            
            switch (ii[3]){
            case 0:
            	t.setComponent(0, 0, 1);
            	t.setComponent(1, 1, 1);
            	t.setComponent(2, 2, 1);
            	break;
            case 1:
            	t.setComponent(0, 0,-1);
            	t.setComponent(1, 1, 1);
            	t.setComponent(2, 2,-1);
            	break;
            case 2:
            	t.setComponent(0, 0,-1);
            	t.setComponent(1, 1,-1);
            	t.setComponent(2, 2,-1);
            	break;
            case 3:
            	t.setComponent(0, 0, 1);
            	t.setComponent(1, 1,-1);
            	t.setComponent(2, 2, 1);
            	break;
            }
      
      	  ((AtomActionTransformed)atomGroupAction.getAtomAction()).setTransformationTensor(t);
          atomGroupAction.actionPerformed(molecule);
 
            atomActionTranslateTo.setDestination((Vector)myLat.site(ii));
            atomActionTranslateTo.actionPerformed(molecule);
            
            if (ii[3] == 0){
            	cellPosition = space.makeVector();
            	cellPosition.E((Vector)myLat.site(ii));
            	
            	
            //remember the coordinate of the cell
            //Loop 8 times over the basis and we can make the cell assignment here!!
            }
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

    /**
     * Returns a lattice with positions the same as those used in the 
     * most recent use of initializeCoordinates.  Includes any scaling
     * or translation applied to fill the space, and thus will not necessarily
     * be the same positions as specified by the lattice given at construction.
     */
    public SpaceLattice getLatticeMemento() {
        return myLat;
    }

    protected final SpaceLattice lattice;
    protected final IndexIteratorSizable indexIterator;
    protected boolean rescalingToFitVolume = true;
    protected final MoleculeActionTranslateTo atomActionTranslateTo;
    protected final MoleculeChildAtomAction atomGroupAction;
    protected MyLattice myLat;
    private static final long serialVersionUID = 2L;
    /**
     * Returns the resizeLatticeToFitVolume flag.
     * 
     * @return boolean
     */
    public boolean isRescalingToFitVolume() {
        return rescalingToFitVolume;
    }

    /**
     * Sets the resizeLatticeToFitVolume flag, which if true indicates that the
     * primitive vectors should be resized to fit the dimensions of the box.
     * Default is true.
     * 
     * @param resizeLatticeToFitVolume
     *            The resizeLatticeToFitVolume to set
     */
    public void setRescalingToFitVolume(boolean resizeLatticeToFitVolume) {
        this.rescalingToFitVolume = resizeLatticeToFitVolume;
    }

    /**
     * Used to store the state of a lattice.
     * 
     * @author nancycribbin, Andrew Schultz, Dr. Kofke
     * 
     */
    public static class MyLattice implements SpaceLattice {

        public MyLattice(SpaceLattice l, Vector latticeScaling, Vector offset) {
            lattice = l;
            this.latticeScaling = latticeScaling;
            this.offset = offset;
            this.site = l.getSpace().makeVector();
        }

        public Space getSpace() {
            return lattice.getSpace();
        }

        public int D() {
            return lattice.D();
        }

        /**
         * Returns the same instance of Vector with each call.
         */
        public Object site(int[] index) {
            site.E((Vector) lattice.site(index));
            site.TE(latticeScaling);
            site.PE(offset);

            return site;
        }

        public double[] getLatticeConstants() {
            double[] lat = lattice.getLatticeConstants();
            for (int i = 0; i < lat.length; i++) {
                lat[i] *= latticeScaling.getX(i);
            }
            return lat;
        }

        final SpaceLattice lattice;
        final public Vector latticeScaling;
        final Vector offset;
        final Vector site;

    }

}
