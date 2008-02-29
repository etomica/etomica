package etomica.spin;

import etomica.api.IVector;
import etomica.atom.AtomSet;
import etomica.atom.IAtomPositioned;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.box.Box;
import etomica.config.Configuration;


/**
 * Configuration used to align all atoms in a spin system, so all
 * point in the same direction.  Spin direction is given by the 
 * atom's position vector.
 *
 * @author David Kofke
 *
 */
public class ConfigurationAligned implements Configuration, java.io.Serializable {

    /**
     * @param space
     */
    public ConfigurationAligned() {
    }

    /**
     * Sets all spins to be aligned in the +x direction
     */
    public void initializeCoordinates(Box box) {
        AtomSet leafAtoms = box.getLeafList();
        for (int i=0; i<leafAtoms.getAtomCount(); i++) {
            IVector spin = ((IAtomPositioned)leafAtoms.getAtom(i)).getPosition();
            spin.E(0.0);
            spin.setX(0,1.0);
        }
    }

    private static final long serialVersionUID = 2L;
}
