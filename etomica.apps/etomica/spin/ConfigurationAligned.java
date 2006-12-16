package etomica.spin;

import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.config.Configuration;
import etomica.phase.Phase;
import etomica.space.Vector;


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
    public void initializeCoordinates(Phase phase) {
        AtomIteratorAllMolecules iterator = new AtomIteratorAllMolecules(phase);
        iterator.reset();
        while(iterator.hasNext()) {
            Vector spin = ((AtomLeaf)iterator.nextAtom()).getCoord().position();
            spin.E(0.0);
            spin.setX(0,1.0);
        }
    }

    private static final long serialVersionUID = 2L;
}
