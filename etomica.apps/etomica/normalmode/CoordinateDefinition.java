package etomica.normalmode;

import etomica.atom.Atom;

/**
 * An interface that defines the real-space generalized coordinates that are
 * summed (over molecules) to form collective coordinates (from which normal
 * coordinates are determined). Typically these generalized coordinates are
 * given by the displacement of the molecule from its nominal lattice position.
 * For non-spherical molecules it will also include elements related to the
 * orientation.
 * 
 * @author Andrew Schultz
 */
public interface CoordinateDefinition {

    /**
     * Returns the number of generalized coordinates associated with each
     * molecule. If no orientational coordinates are involved, this value is
     * typically equal to the space dimension.
     */
    public int getCoordinateDim();

    /**
     * Notifies the coord definition how many atoms will be tracked. The
     * |index| parameter in other methods must not exceed |numAtoms-1|.
     */
    public void setNumAtoms(int numAtoms);

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
    public void calcU(Atom molecule, int index, double[] u);

    /**
     * Initializes the CoordinateDefinition for the given molecule and
     * associates the molecule with the given index. Typically this will be
     * called when the molecule is in a nominal position and orientation, and
     * the generalized coordinates for the molecule will be defined with respect
     * to this nominal case.
     */
    public void initNominalU(Atom molecule, int index);

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
    public void setToU(Atom molecule, int index, double[] u);
}
