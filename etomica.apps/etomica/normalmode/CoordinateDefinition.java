package etomica.normalmode;

import java.io.Serializable;

import etomica.api.IAtomSet;
import etomica.api.IAtom;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IVector;

import etomica.action.AtomActionTranslateTo;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomTypeMolecule;
import etomica.atom.AtomsetArrayList;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.config.Conformation;
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

    public CoordinateDefinition(IBox box, int coordinateDim, Primitive primitive, Space _space) {
        this(box, coordinateDim, primitive, new BasisMonatomic(_space), _space);
    }
    
    public CoordinateDefinition(IBox box, int coordinateDim, Primitive primitive, Basis basis, Space _space) {
        this.coordinateDim = coordinateDim;
        this.primitive = primitive;
        this.basis = basis;
        this.box = box;
        this.space = _space;
        lattice = new BravaisLatticeCrystal(primitive, basis);

        atomActionTranslateTo = new AtomActionTranslateTo(lattice.getSpace());

        
    }
    
    public void initializeCoordinates(int[] nCells) {
        AtomIteratorAllMolecules atomIterator = new AtomIteratorAllMolecules(box);
        IAtomSet moleculeList = box.getMoleculeList();
        if (moleculeList.getAtomCount() == 0) {
            throw new RuntimeException("There are no atoms yet!");
        }

        int basisSize = lattice.getBasis().getScaledCoordinates().length;

        IVector offset = lattice.getSpace().makeVector();
        IVector[] primitiveVectors = primitive.vectors();
        for (int i=0; i<primitiveVectors.length; i++) {
            offset.PEa1Tv1(nCells[i],primitiveVectors[i]);
        }
        offset.TE(-0.5);
        
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
        IVector position = lattice.getSpace().makeVector();
        AtomArrayList currentList = null;
        for (int iMolecule = 0; iMolecule<moleculeList.getAtomCount(); iMolecule++) {
            IMolecule molecule = (IMolecule)moleculeList.getAtom(iMolecule);
            // initialize coordinates of child atoms
            Conformation config = ((AtomTypeMolecule)molecule.getType()).getConformation();
            config.initializePositions(molecule.getChildList());

            int[] ii = indexIterator.next();
            position.E((IVector)lattice.site(ii));
            position.PE(offset);
            
            atomActionTranslateTo.setDestination(position);
            atomActionTranslateTo.actionPerformed(molecule);

            if (ii[space.D()] == 0) {
                if (iCell > -1) {
                    initNominalU(cells[iCell].molecules);
                }
                // new cell
                iCell++;
                currentList = new AtomArrayList(basisSize);
                cells[iCell] = new BasisCell(new AtomsetArrayList(currentList), lattice.getSpace().makeVector());
                cells[iCell].cellPosition.E(position);
            }
            currentList.add(molecule);
        }
        
        initNominalU(cells[totalCells-1].molecules);
        
        siteManager = new AtomAgentManager(new SiteSource(space), box);
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
     * @param molecule
     *            The molecule of interest, which should be those forming a unit cell of the lattice
     */
    protected abstract double[] calcU(IAtomSet molecules);

    /**
     * Initializes the CoordinateDefinition for the given molecule and
     * associates the molecule with the given index. Typically this will be
     * called when the molecule is in a nominal position and orientation, and
     * the generalized coordinates for the molecule will be defined with respect
     * to this nominal case.
     */
    protected abstract void initNominalU(IAtomSet molecules);

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
    public abstract void setToU(IAtomSet molecules, double[] newU);

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
    public void calcT(IVector k, double[] realT, double[] imaginaryT) {
        for (int i = 0; i < coordinateDim; i++) {
            realT[i] = 0;
            imaginaryT[i] = 0;
        }

        // sum T over atoms
        for (int iCell = 0; iCell<cells.length; iCell++) {
            BasisCell cell = cells[iCell];
            IAtomSet molecules = cell.molecules;
            double[] u = calcU(molecules);
            IVector latticePosition = cell.cellPosition;
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

    public IBox getBox() {
        return box;
    }

    public IVector getLatticePosition(IAtom atom) {
        return (IVector)siteManager.getAgent(atom);
    }
    
    public BasisCell[] getBasisCells() {
        return cells;
    }

    protected final int coordinateDim;
    protected final IBox box;
    protected AtomAgentManager siteManager;
    protected final BravaisLatticeCrystal lattice;
    protected final Primitive primitive;
    protected final Basis basis;
    protected final AtomActionTranslateTo atomActionTranslateTo;
    protected BasisCell[] cells;
    protected final Space space;
    
    protected static class SiteSource implements AgentSource, Serializable {
        
        public SiteSource(Space space) {
            this.space = space;
        }
        public Class getAgentClass() {
            return IVector.class;
        }
        public Object makeAgent(IAtom atom) {
            IVector vector = space.makeVector();
            vector.E(atom.getType().getPositionDefinition().position(atom));
            return vector;
        }
        public void releaseAgent(Object agent, IAtom atom) {
            //nothing to do
        }

        private final Space space;
        private static final long serialVersionUID = 1L;
    }
    
    public static class BasisCell implements Serializable {
        public BasisCell(IAtomSet molecules, IVector cellPosition) {
            this.molecules = molecules;
            this.cellPosition = cellPosition;
        }
        
        private static final long serialVersionUID = 1L;
        public final IAtomSet molecules;
        public final IVector cellPosition;
    }

}
