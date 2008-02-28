package etomica.models.water;

import etomica.atom.DipoleSource;
import etomica.atom.IMolecule;
import etomica.atom.MoleculeOriented;
import etomica.space.IVector;
import etomica.space.Space;

/**
 * Implementation of DipoleSource that can handle water molecules that are
 * MoleculeOriented.  The orientation is taken from the molecule rather than
 * calculated from scratch.
 * 
 * The dipole is assumed to point in the "secondary" direction held by the
 * orientation.  This is the correct direction for AtomWater3P.  The magnitude
 * of the dipole must be given to this class (via setDipoleStrength since that
 * depends on the potential in use.
 *
 * @author Andrew Schultz
 */
public class DipoleSourceWater implements DipoleSource {

    public DipoleSourceWater(Space space) {
        dipole = space.makeVector();
    }

    /**
     * Sets the strength (magnitude) of the dipole.
     */
    public void setDipoleStrength(double newDipoleStrength) {
        dipoleStrength = newDipoleStrength;
    }

    /**
     * Returns the dipole strength.
     */
    public double getDipoleStrength() {
        return dipoleStrength;
    }
    
    public IVector getDipole(IMolecule molecule) {
        // assume dipole points in the secondary orientation direction
        MoleculeOriented orientedMolecule = (MoleculeOriented)molecule;
        dipole.E(orientedMolecule.getOrientation().getSecondaryDirection());
        dipole.TE(dipoleStrength);
        return dipole;
    }

    protected double dipoleStrength;
    protected final IVector dipole;
}
