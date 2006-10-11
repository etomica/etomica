package etomica.models.hexane;

import etomica.atom.Atom;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.lattice.Primitive;
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
        // instance instead of the actual boundary so that boundary changes
        // don't affect us (since we work off original coordinates)
        bdry = new BoundaryRectangularPeriodic(phase.space(), 1);
        bdry.setDimensions(phase.getBoundary().getDimensions());

        dim = phase.space().D();
        maxes = new int[dim];
        mins = new int[dim];
        jumpCount = new int[dim];
        halfMax = new int[dim];
        temp = phase.space().makeVector();

        inverter = phase.space().makeTensor();
        inverter.E(prim.vectors());
        inverter.inverse();
        
        indicesV = new Vector[phase.getSpeciesMaster().getMaxGlobalIndex() + 1];

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
        r0.E(firstatom.type.getPositionDefinition().position(firstatom));

        temp.E(0.0);
        while (aim.hasNext()) {
            Atom molecule = aim.nextAtom();
            temp.E(molecule.type.getPositionDefinition().position(molecule));

            indicesV[molecule.getGlobalIndex()] = (Vector) temp.clone();
            temp.ME(r0);

            bdry.nearestImage(temp);

            tempI = calculateTheseIndices(temp);

            for (int i = 0; i < dim; i++) {
                if (maxes[i] < tempI[i]) {
                    maxes[i] = tempI[i];
                }
                if (mins[i] > tempI[i]) {
                    mins[i] = tempI[i];
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
            if (iJump[i] % 2 == 0) {
                // if jump is even, -half and +half are equivalent and we'll
                // need to detect that below and fix it.
                halfMax[i] = iJump[i] / 2;
            }
            else {
                // we don't actually need halfMax for this case
                // set it to something bogus
                halfMax[i] = iJump[i]*2;
            }
        }
        maxLength = jumpCount[0] * iJump[0];
    }

    /**
     * Returns the Miller indices of the vector argument
     */
    // nan is this okay after the atom has moved?
    private int[] calculateTheseIndices(Vector v) {
        v.transform(inverter);

        int[] w = new int[v.D()];
        for (int i = 0; i < v.D(); i++) {
            w[i] = (int) Math.round(v.x(i));
        }

        return w;
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
        temp.Ev1Mv2(indicesV[atom1.getGlobalIndex()],
                indicesV[atom0.getGlobalIndex()]);

        bdry.nearestImage(temp);

        tempI = calculateTheseIndices(temp);

        // Now we make sure the leading Miller index is either 0 or positive,
        // which
        // is a definite framework for how the vectors point.
        orderMillerIndices(tempI);

        //XXX ick.  If -half and +half are identical, then make -half become
        //+half
        for (int i=0; i<dim; i++) {
            if (tempI[i] == -halfMax[i]) {
                tempI[i] = halfMax[i];
            }
        }

        // Calculate the bin number based on the Miller indices and the atom
        // numbers in the molecule.
        return calculateBinNumber(tempI);
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
            if (m[i] == halfMax[i] || m[i] == -halfMax[i]) {
                // on the edge.  this can get flipped later if it's negative
                // after we do our thing.
                continue;
            }
            if (m[i] > 0) {
                //all is well
                break;
            }
            else if (m[i] < 0) {
                //all is backwards.  fix it
                flipflag = true;
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

    public boolean isFlipFlag() {
        return flipflag;
    }

    Phase phase; // Contains the atom information.

    int dim; // The dimensions of the system

    private int[] maxes; // The maximum values in each physical direction
    private int[] mins;
    private int[] halfMax;
    private int[] jumpCount;

    private int maxLength; // The maximum number of index values possible. May

    // be needed outside of this class.

    Primitive prim; // The primitives used to define the space

    Vector temp; // Temporary storage space.

    Vector[] tempV; // Temporary storage space.

    int[] tempI; // Temporary storage space.

    Vector[] indicesV; // The position-in-space vector between the firstatom

    // and the atom.

    Tensor inverter;

    Boundary bdry;

    boolean flipflag;
}