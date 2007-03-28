package etomica.normalmode;

import etomica.atom.Atom;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.phase.Phase;
import etomica.space.IVector;

/**
 * An interface that defines the real-space generalized coordinates that are
 * summed (over molecules) to form collective coordinates (from which normal
 * coordinates are determined). Typically these generalized coordinates are
 * given by the displacement of the molecule from its nominal lattice position.
 * For non-spherical molecules it will also include elements related to the
 * orientation.
 * 
 * @author Andrew Schultz, David Kofke
 */
public abstract class CoordinateDefinition {

    public CoordinateDefinition(int coordinateDim) {
        this.coordinateDim = coordinateDim;
        u = new double[coordinateDim];
        iterator = new AtomIteratorAllMolecules();
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
     * Calculates the generalized coordinates for the given molecule in its
     * current position and orientation.
     * 
     * @param molecule
     *            The molecule of interest
     * @param index
     *            The index for the molecule as specified via initNominalU
     * @param u
     *            Upon return, the atom's generalized coordinates. |u| must be
     *            of length getCoordinateDim()
     */
    public abstract void calcU(Atom molecule, int index, double[] u);

    /**
     * Initializes the CoordinateDefinition for the given molecule and
     * associates the molecule with the given index. Typically this will be
     * called when the molecule is in a nominal position and orientation, and
     * the generalized coordinates for the molecule will be defined with respect
     * to this nominal case.
     */
    public abstract void initNominalU(Atom molecule, int index);

    /**
     * Set the molecule to a position and orientation that corresponds to the
     * given generalized coordinate. |u| must be of length getNormalDim()
     * 
     * @param molecule
     *            The molecule of interest
     * @param index
     *            The index for the molecule as specified via initNominalU
     * @param u
     *            The generalized coordinate that defines the position and
     *            orientation to which the molecule will be set by this method.
     */
    public abstract void setToU(Atom molecule, int index, double[] u);

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
    public void calcT(IVector k, double[] realT, double[] imaginaryT) {
        for (int i = 0; i < coordinateDim; i++) {
            realT[i] = 0;
            imaginaryT[i] = 0;
        }
        iterator.reset();
        int index = 0;
        // sum T over atoms
        while (iterator.hasNext()) {
            Atom atom = iterator.nextAtom();
            calcU(atom, index, u);
            double kR = k.dot(latticePositions[index]);
            double coskR = Math.cos(kR);
            double sinkR = Math.sin(kR);
            for (int i = 0; i < coordinateDim; i++) {
                realT[i] += coskR * u[i];
                imaginaryT[i] += sinkR * u[i];
            }
            index++;
        }

        for (int i = 0; i < coordinateDim; i++) {
            realT[i] /= sqrtN;
            imaginaryT[i] /= sqrtN;
        }

    }

    public Phase getPhase() {
        return phase;
    }

    public void setPhase(Phase phase) {
        this.phase = phase;
        N = phase.getSpeciesMaster().moleculeCount();
        sqrtN = Math.sqrt((double) N);
        latticePositions = new IVector[N];

        iterator.setPhase(phase);
        iterator.reset();
        int index = 0;
        while (iterator.hasNext()) {
            latticePositions[index] = phase.getSpace().makeVector();
            Atom atom = iterator.nextAtom();
            latticePositions[index].E(atom.getType().getPositionDefinition()
                    .position(atom));
            index++;
        }

        // record nominal position of each atom
        nominalU = new double[iterator.size()][getCoordinateDim()];
        index = 0;
        iterator.reset();
        while (iterator.hasNext()) {
            initNominalU(iterator.nextAtom(), index++);
        }

    }

    public IVector[] getLatticePositions() {
        return latticePositions;
    }

    /**
     * Notifies the coord definition how many molecules will be tracked. The |index|
     * parameter in other methods must not exceed |numAtoms-1|.
     */
    public abstract void setNumAtoms(int numAtoms);

    protected final int coordinateDim;

    protected double[][] nominalU;

    private final double[] u;

    protected Phase phase;

    private int N;

    private double sqrtN;

    private IVector[] latticePositions;

    private final AtomIteratorAllMolecules iterator;

}
