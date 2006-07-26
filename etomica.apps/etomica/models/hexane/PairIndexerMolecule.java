package etomica.models.hexane;

import etomica.atom.Atom;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomPair;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.atom.iterator.AtomIteratorTree;
import etomica.lattice.Primitive;
import etomica.phase.Phase;
import etomica.space.Boundary;
import etomica.space.Tensor;
import etomica.space.Vector;

/**
 * Given molecules in the form of an AtomPair, this class returns an integer index 
 * for it, based on the positions of the molecules.
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
        this.bdry = (Boundary) phase.getBoundary().clone();

        maxes = new int[dim];
        dim = phase.space().D();
        temp = phase.space().makeVector();

        inverter = phase.space().makeTensor();
        inverter.E(prim.vectors());
        inverter.inverse();

        calculateAllAtomIndices();
        // Make the maximum length of the storage array
        maxLength = 1;
        for (int i = 0; i < dim; i++) {
            maxLength *= getMaxes(i);
        }
    }

    /**
     * Assign the molecule index to each individual atom.
     */
    private void calculateAllAtomIndices() {
        indices = new int[phase.getSpeciesMaster().getMaxGlobalIndex() + 1][];

        AtomIteratorAllMolecules aim = new AtomIteratorAllMolecules(phase);
        aim.reset();
        Atom firstatom = (Atom) aim.peek();
        Vector r0 = phase.space().makeVector();
        r0.E(firstatom.type.getPositionDefinition().position(firstatom));

        temp.E(0.0);
        AtomIteratorTree ait = new AtomIteratorTree();
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
            }

        }
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
    public int getBin(AtomPair atompair) {

        // All these calculations are done in "primitive units"

        // create the vector between the atoms.
        temp.Ev1Mv2(indicesV[atompair.getAtom(1).getGlobalIndex()],
                indicesV[atompair.getAtom(0).getGlobalIndex()]);

        bdry.nearestImage(temp);

        tempI = calculateTheseIndices(temp);

        // Now we make sure the leading Miller index is either 0 or positive,
        // which
        // is a definite framework for how the vectors point.
        orderMillerIndices(tempI);

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
        // Instead of a lot of confusing loops, why not just use a switch
        // statment?
        switch (dim) {
        case 1: {
            if (m[0] < 0) {
                m[0] *= -1;
                flipflag = true;
            }
        }
            break;
        case 2: {
            if (m[0] < 0) {
                m[0] *= -1;
                m[1] *= -1;
                flipflag = true;
            } else if (m[0] == 0) {
                if (m[1] < 0) {
                    m[1] *= -1;
                    flipflag = true;
                }
            }
        }
            break;
        case 3: {
            if (m[0] < 0) {
                m[0] *= -1;
                m[1] *= -1;
                m[2] *= -1;
                flipflag = true;
            } else if (m[0] == 0) {
                if (m[1] < 0) {
                    m[1] *= -1;
                    m[2] *= -1;
                    flipflag = true;
                } else if (m[1] == 0) {
                    if (m[2] < 0) {
                        m[2] *= -1;
                        flipflag = true;
                    }
                }
            }
        }
            break;
        }// end switch statement

    }// end method

    /**
     * Given the vector between the original locations of the atoms, in
     * primitive units, and the atoms themselves, calculates the bin number that
     * information should be stored in.
     */
    private int calculateBinNumber(int[] tin) {
        int t = 0;

        switch (dim) {
        case 1: {
            t = tin[0];
        }
            break;
        case 2: {
            t = tin[0] * getMaxes(1) + tin[1];

        }
            break;
        case 3: {
            t = tin[0] * (getMaxes(1) + getMaxes(2)) + tin[1] * getMaxes(2)
                    + tin[2];
        }
            break;
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

    private int maxLength; // The maximum number of index values possible. May
                            // be needed outside of this class.

    Primitive prim; // The primitives used to define the space

    Vector temp; // Temporary storage space.

    Vector[] tempV; // Temporary storage space.

    int[] tempI; // Temporary storage space.

    int[][] indices; // The first index is the global index of the atom.

    Vector[] indicesV;  // The position-in-space vector between the firstatom
                        // and the atom.

    Tensor inverter;

    Boundary bdry;

    boolean flipflag;
}