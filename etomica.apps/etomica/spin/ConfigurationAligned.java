package etomica.spin;

import etomica.atom.IAtomPositioned;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.config.Configuration;
import etomica.box.Box;
import etomica.space.IVector;


/**
 * Configuration used to align all atoms in a spin system, so all
 * point in the same direction.  Spin direction is given by the 
 * atom's position vector.
 *
 * @author David Kofke
 *
 */
public class ConfigurationAligned extends Configuration {

    /**
     * @param space
     */
    public ConfigurationAligned() {
    }

    /**
     * Sets all spins to be aligned in the +x direction
     */
    public void initializeCoordinates(Box box) {
        AtomIteratorAllMolecules iterator = new AtomIteratorAllMolecules(box);
        iterator.reset();
        for (IAtomPositioned atom = (IAtomPositioned)iterator.nextAtom(); atom != null;
             atom = (IAtomPositioned)iterator.nextAtom()) {
            IVector spin = atom.getPosition();
            spin.E(0.0);
            spin.setX(0,1.0);
        }
    }

    private static final long serialVersionUID = 2L;
}
