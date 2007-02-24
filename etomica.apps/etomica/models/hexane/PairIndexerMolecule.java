package etomica.models.hexane;

import etomica.atom.Atom;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.lattice.crystal.Primitive;
import etomica.phase.Phase;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Tensor;
import etomica.space.Vector;

/**
 * Given molecules in the form of an AtomPair, this class returns an integer
 * index for it, based on the positions of the molecules.
 * 
 * 
 * MOLECULE LEVEL- DOES NOT DEAL WITH ATOMS
 * 
 * @author nancycribbin
 */

// nan MOLDECULES ONLY!!!!!!!!!!!!!
public class PairIndexerMolecule {
    // make sure you aren't losing the redundancy:
    // fcc lattice = simple cell with basis of 4 atoms.
    // says: I've got a primitive and phase.
    // I need to look at the molecule layer.
    // AtomPositionDefinition comes from AtomFactory of 'cule which is held by
    // the type or the species
    //

    public PairIndexerMolecule(Phase ph, Primitive pr) {
        // Initialize everything we can.
        this.phase = ph;
        this.prim = pr.copy();
        // clone because we don't know how the boundary will change, and we want
        // an undisturbed
        // original.
//        this.bdry = (Boundary) phase.getBoundary().clone();
        
        // assume we have a rectangular box.  make our own little boundary
        // instance instead of using actual boundary so that boundary changes
        // don't affect us (since we work off original coordinates)
        bdry = new BoundaryRectangularPeriodic(phase.space(), 1);
        bdry.setDimensions(phase.getBoundary().getDimensions());

        dim = phase.space().D();
        maxes = new int[dim];
        mins = new int[dim];
        jumpCount = new int[dim];
        temp = phase.space().makeVector();
        indices = new int[dim];

        inverter = phase.space().makeTensor();
        inverter.E(prim.vectors());
        inverter.inverse();
        
        latticeSites = new Vector[phase.getSpeciesMaster().getMaxGlobalIndex() + 1];

        calculateAllAtomIndices();
    }

    /**
     * Assign the molecule index to each individual atom.
     */
    private void calculateAllAtomIndices() {
        for (int i=0; i<dim; i++) {
            maxes[i] = Integer.MIN_VALUE;
            mins[i] = Integer.MAX_VALUE;
        }

        AtomIteratorAllMolecules aim = new AtomIteratorAllMolecules(phase);
        aim.reset();
        Atom firstatom = (Atom) aim.peek();
        Vector r0 = phase.space().makeVector();
        r0.E(firstatom.getType().getPositionDefinition().position(firstatom));

        temp.E(0.0);
        while (aim.hasNext()) {
            Atom molecule = aim.nextAtom();
            temp.E(molecule.getType().getPositionDefinition().position(molecule));

            latticeSites[molecule.getGlobalIndex()] = phase.space().makeVector();
            latticeSites[molecule.getGlobalIndex()].E(temp);
            temp.ME(r0);

            bdry.nearestImage(temp);
            
            // flip the vector to ensure symmetry
            //XXX we can only do this with spheres!!!
            flipVector(temp);

            calculateTheseIndices(temp);
            orderMillerIndices(indices);

            for (int i = 0; i < dim; i++) {
                if (maxes[i] < indices[i]) {
                    maxes[i] = indices[i];
                }
                if (mins[i] > indices[i]) {
                    mins[i] = indices[i];
                }
            }
        }

        int[] iJump = new int[dim];
        for (int i=dim-1; i>-1; i--) {
            // iJump is how many indices exist for dimension i
            iJump[i] = maxes[i] - mins[i] + 1;
            if (i==dim-1) {
                // jumpCount is how many bins we increase each time index i 
                // increases by 1 with all other indices held fixed
                jumpCount[i] = 1;
            }
            else {
                jumpCount[i] = jumpCount[i+1] * iJump[i+1];
            }
        }
        maxLength = jumpCount[0] * iJump[0];
    }
    
    private void flipVector(Vector dr) {
        int shouldBeFlipped = 0; // 1=yes, -1=no
        for (int i=0; i<dr.D(); i++) {
            if (Math.abs(dr.x(i) - 0.5*bdry.getDimensions().x(i)) < tol || 
                Math.abs(dr.x(i) + 0.5*bdry.getDimensions().x(i)) < tol) {
                // we're on the edge, put ourselves on the right edge
                dr.setX(i, 0.5*bdry.getDimensions().x(i));
            }
            else if (shouldBeFlipped == -1) {
                // we already hit a positive, so don't flip this
            }
            else if (shouldBeFlipped == 1) {
                // we already hit a negative, so flip this
                dr.setX(i, -dr.x(i));
            }
            else if (Math.abs(dr.x(i)) < tol) {
                // this is a 0, do nothing
            }
            else if (dr.x(i) < 0) {
                // this is a negative, flip it and all the rest
                shouldBeFlipped = 1;
                dr.setX(i, -dr.x(i));
            }
            else {
                // this is a positive, don't flip it or the rest
                shouldBeFlipped = -1;
            }
        }
    }

    /**
     * Calculates the Miller indices of the vector argument
     */
    // nan is this okay after the atom has moved?
    private void calculateTheseIndices(Vector v) {
        // transformation is done in place
        v.transform(inverter);

        for (int i = 0; i < v.D(); i++) {
            indices[i] = (int) Math.round(v.x(i));
        }
    }

    protected int[] getMaxes() {
        return maxes;
    }

    protected int getMaxes(int i) {
        return maxes[i];
    }

    /**
     * Returns the appropriate bin number for storing information for a given
     * pair of atoms
     */
    public int getBin(Atom atom0, Atom atom1) {

        // All these calculations are done in "primitive units"

        // create the vector between the atoms.
        temp.Ev1Mv2(latticeSites[atom1.getGlobalIndex()],
                latticeSites[atom0.getGlobalIndex()]);

        bdry.nearestImage(temp);
        flipVector(temp);

        calculateTheseIndices(temp);

        // Now we make sure the leading Miller index is either 0 or positive,
        // which
        // is a definite framework for how the vectors point.
        orderMillerIndices(indices);

        for (int i=0; i<dim; i++) {
            if (indices[i] < mins[i] || indices[i] > maxes[i]) {
                System.out.println("indices="+indices[0]+" "+indices[1]+" "+indices[2]);
                throw new RuntimeException("index["+i+"]="+indices[i]+" is out of bounds "+temp);
            }
        }

        // Calculate the bin number based on the Miller indices and the atom
        // numbers in the molecule.
        return calculateBinNumber(indices);
    }

    /**
     * This method makes sure that the first non-zero Miller index is positive.
     * It is a way to make sure all the vectors between 2 lattice points are
     * oriented in the same direction, so that the bin which is calculated will
     * be correct.
     */
    // Note that this method operates on and may change its argument.
    private void orderMillerIndices(int[] m) {
        // find the first non-zero index.
        for (int i = 0; i < dim; i++) {
            if (m[i] > 0) {
                //all is well
                break;
            }
            else if (m[i] < 0) {
                //all is backwards.  fix it
                for (int j = i; j < dim; j++) {
                    m[j] = -m[j];
                }
                break;
            }
            // else it's 0, so go on to the next one
        }
    }// end method

    /**
     * Given the vector between the original locations of the atoms, in
     * primitive units, and the atoms themselves, calculates the bin number that
     * information should be stored in.
     */
    private int calculateBinNumber(int[] tin) {
        int t = 0;
        for (int i=0; i<dim; i++) {
            t += tin[i] * jumpCount[i];
        }
        if (t < 0) {
            throw new RuntimeException("bin can't be negative "+t+" "+tin[0]+" "+tin[1]+" "+tin[2]);
        }
        return t;
    }

    public int getMaxLength() {
        return maxLength;
    }

    /**
     * Returns the array of Vectors corresponding to the original position of 
     * each Atom relative to the first Atom.  These should correspond to the 
     * actual lattice sites the Atoms were placed at.
     */
    public Vector[] getLatticeSites() {
        return latticeSites;
    }

    /**
     * Contains the atom information.
     */
    private Phase phase;

    /**
     * The dimensions of the system
     */
    private int dim;

    /**
     * The maximum values in each physical direction
     */
    private int[] maxes;
    private int[] mins;
    private int[] jumpCount;

    /**
     * The maximum number of index values possible.
     */
    // Kept as a field in case it's needed outside this class
    private int maxLength;
    
    /**
     * The primitives used to define the lattice
     */
    private final Primitive prim; 

    /**
     * Temporary storage space for a vector
     */
    private final Vector temp;

    /**
     * Storage space to put the indices
     */
    private final int[] indices;

    /**
     * The original position-in-space vector of each atom
     */
    private final Vector[] latticeSites; 

    private final Tensor inverter;

    private final Boundary bdry;
    private double tol = 1.e-5;
}