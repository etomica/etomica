/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import java.io.Serializable;

import etomica.action.MoleculeActionTranslateTo;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.atom.IMolecule;
import etomica.atom.IMoleculeList;
import etomica.space.Vector;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.MoleculeArrayList;
import etomica.atom.MoleculeListWrapper;
import etomica.atom.iterator.MoleculeIteratorAllMolecules;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.IndexIteratorRectangular;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisMonatomic;
import etomica.lattice.crystal.Primitive;
import etomica.space.Space;

/**
 * An abstract class that defines the real-space generalized coordinates that are
 * summed (over molecules) to form collective coordinates (from which normal
 * coordinates are determined). Typically these generalized coordinates are
 * given by the displacement of the molecule from its nominal lattice position.
 * For non-spherical molecules it will also include elements related to the
 * orientation.  This class provides methods to compute the k-dependent collective generalized
 * coordinates, obtained by Fourier-summing the generalized coordinates over the atoms
 * in a box.
 * 
 * @author Andrew Schultz, David Kofke
 */
public abstract class CoordinateDefinition {

    public CoordinateDefinition(Box box, int coordinateDim, Primitive primitive, Space _space) {
        this(box, coordinateDim, primitive, new BasisMonatomic(_space), _space);
    }
    
    public CoordinateDefinition(Box box, int coordinateDim, Primitive primitive, Basis basis, Space _space) {
        this.coordinateDim = coordinateDim;
        this.primitive = primitive;
        this.basis = basis;
        this.box = box;
        this.space = _space;
        lattice = new BravaisLatticeCrystal(primitive, basis);

        atomActionTranslateTo = new MoleculeActionTranslateTo(lattice.getSpace());
    }
    
    public void initializeCoordinates(int[] nCells) {
        MoleculeIteratorAllMolecules atomIterator = new MoleculeIteratorAllMolecules(box);
        IMoleculeList moleculeList = box.getMoleculeList();
        if (moleculeList.getMoleculeCount() == 0) {
            throw new RuntimeException("There are no atoms yet!");
        }

        int basisSize = lattice.getBasis().getScaledCoordinates().length;
        Vector minBasis = space.makeVector();
        minBasis.E(1e10);
        Vector maxBasis = space.makeVector();
        maxBasis.E(-1e10);
        for (int i=0; i<basisSize; i++) {
            Vector p = lattice.getBasis().getScaledCoordinates()[i];
            for (int j=0; j<minBasis.getD(); j++) {
                if (minBasis.getX(j) > p.getX(j)) {
                    minBasis.setX(j, p.getX(j));
                }
                if (maxBasis.getX(j) < p.getX(j)) {
                    maxBasis.setX(j, p.getX(j));
                }
            }
        }
        Vector basisOffset = space.makeVector();
        basisOffset.E(maxBasis);
        basisOffset.PE(minBasis);
        basisOffset.TE(-0.5);
        basisOffset.PE(0.5);

        Vector offset = lattice.getSpace().makeVector();
        Vector[] primitiveVectors = primitive.vectors();
        for (int i=0; i<primitiveVectors.length; i++) {
            offset.PEa1Tv1(-0.5*nCells[i],primitiveVectors[i]);
            offset.PEa1Tv1(basisOffset.getX(i), primitiveVectors[i]);
        }
        
        IndexIteratorRectangular indexIterator = new IndexIteratorRectangular(space.D()+1);
        int[] iteratorDimensions = new int[space.D()+1];
        
        
        System.arraycopy(nCells, 0, iteratorDimensions, 0, nCells.length);
        iteratorDimensions[nCells.length] = basisSize;
        indexIterator.setSize(iteratorDimensions);

        int totalCells = 1;
        for (int i=0; i<nCells.length; i++) {
            totalCells *= nCells[i];
        }
        
        cells = new BasisCell[totalCells];
        int iCell = -1;
        // Place molecules
        atomIterator.reset();
        indexIterator.reset();
        Vector position = lattice.getSpace().makeVector();
        MoleculeArrayList currentList = null;
        for (int iMolecule = 0; iMolecule<moleculeList.getMoleculeCount(); iMolecule++) {
            IMolecule molecule = moleculeList.getMolecule(iMolecule);
            // initialize coordinates of child atoms
            molecule.getType().initializeConformation(molecule);

            int[] ii = indexIterator.next();
            position.E((Vector)lattice.site(ii));
            position.PE(offset);
            
            atomActionTranslateTo.setDestination(position);
            atomActionTranslateTo.actionPerformed(molecule);

            if (ii[space.D()] == 0) {
                if (iCell > -1) {
                    initNominalU(cells[iCell].molecules);
                }
                // new cell
                iCell++;
                currentList = new MoleculeArrayList(basisSize);
                cells[iCell] = new BasisCell(new MoleculeListWrapper(currentList), lattice.getSpace().makeVector());
                cells[iCell].cellPosition.E(position);
            }
            currentList.add(molecule);
        }
        
        initNominalU(cells[totalCells-1].molecules);
        
        siteManager = new AtomLeafAgentManager<Vector>(new SiteSource(space), box, Vector.class);
    }

    

    /**
     * Returns the number of generalized coordinates associated with each
     * molecule. If no orientational coordinates are involved, this value is
     * typically equal to the space dimension.
     */
    public int getCoordinateDim() {
        return coordinateDim;
    }

    /**
     * Calculates the generalized coordinates for the given molecules in their
     * current position and orientation.  The result is stored in the |u| field.
     * 
     * @param molecules
     *            The molecules of interest, which should be those forming a unit cell of the lattice
     */
    public abstract double[] calcU(IMoleculeList molecules);

    /**
     * Initializes the CoordinateDefinition for the given molecule and
     * associates the molecule with the given index. Typically this will be
     * called when the molecule is in a nominal position and orientation, and
     * the generalized coordinates for the molecule will be defined with respect
     * to this nominal case.
     */
    protected abstract void initNominalU(IMoleculeList molecules);

    /**
     * Set all the molecules in a cell to a position and orientation that corresponds to the
     * given generalized coordinate. |u| must be of length getCoordinateDim()
     * 
     * @param molecules
     *            The molecules of interest
     * @param newU
     *            The generalized coordinate that defines the position and
     *            orientation to which the molecules will be set by this method.
     */
    public abstract void setToU(IMoleculeList molecules, double[] newU);

    /**
     * Calculates the complex "T vector", which is collective coordinate given
     * by the Fourier sum (over atoms) of the generalized coordinate vector.
     * 
     * @param k
     *            the wave vector
     * @param realT
     *            outputs the real component of the T vector
     * @param imaginaryT
     *            outputs the imaginary component of the T vector
     */
    //in principle this should be returning Complex[] and not returning the values through the args
    public void calcT(Vector k, double[] realT, double[] imaginaryT) {
        for (int i = 0; i < coordinateDim; i++) {
            realT[i] = 0;
            imaginaryT[i] = 0;
        }
        
        // sum T over atoms
        for (int iCell = 0; iCell<cells.length; iCell++) {
        	
            BasisCell cell = cells[iCell];
            IMoleculeList molecules = cell.molecules;
            double[] u = calcU(molecules);
            Vector latticePosition = cell.cellPosition;
            double kR = k.dot(latticePosition);
            double coskR = Math.cos(kR);
            double sinkR = Math.sin(kR);
            for (int i = 0; i < coordinateDim; i++) {
                realT[i] += coskR * u[i];
                imaginaryT[i] -= sinkR * u[i];
            }
        }
        
        double sqrtCells = Math.sqrt(cells.length);
        for (int i = 0; i < coordinateDim; i++) {
            realT[i] /= sqrtCells;
            imaginaryT[i] /= sqrtCells;
        }

    }

    public Box getBox() {
        return box;
    }

    public Vector getLatticePosition(IAtom atom) {
        // this impl only handles leaf atoms.  subclasses might override this
        // method and handle IMolecules.
        return siteManager.getAgent(atom);
    }
    
    public BasisCell[] getBasisCells() {
        return cells;
    }

    public Primitive getPrimitive() {
        return primitive;
    }

    public Basis getBasis() {
        return basis;
    }

    public AtomLeafAgentManager<Vector> getSiteManager() {
        return siteManager;
    }

    protected final int coordinateDim;
    protected final Box box;
    protected AtomLeafAgentManager<Vector> siteManager;
    protected final BravaisLatticeCrystal lattice;
    protected final Primitive primitive;
    protected final Basis basis;
    protected final MoleculeActionTranslateTo atomActionTranslateTo;
    protected BasisCell[] cells;
    protected final Space space;
    
    protected static class SiteSource implements AtomLeafAgentManager.AgentSource<Vector> {
        
        public SiteSource(Space space) {
            this.space = space;
        }
        public Vector makeAgent(IAtom atom, Box agentBox) {
            Vector vector = space.makeVector();
            vector.E(atom.getPosition());
            return vector;
        }
        public void releaseAgent(Vector agent, IAtom atom, Box agentBox) {
            //nothing to do
        }

        private final Space space;
    }
    
    public static class BasisCell implements Serializable {
        public BasisCell(IMoleculeList molecules, Vector cellPosition) {
            this.molecules = molecules;
            this.cellPosition = cellPosition;
        }
        
        private static final long serialVersionUID = 1L;
        public final IMoleculeList molecules;
        public final Vector cellPosition;
    }

}
