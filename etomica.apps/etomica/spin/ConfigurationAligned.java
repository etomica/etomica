package etomica.spin;

import etomica.atom.AtomArrayList;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.config.Configuration;
import etomica.space.Space;
import etomica.space.Vector;


/**
 * Configuration used to align all atoms in a spin system, so all
 * point in the same direction.  Spin direction is given by the 
 * atom's position vector.
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on May 22, 2005 by kofke
 */
public class ConfigurationAligned extends Configuration {

    /**
     * @param space
     */
    public ConfigurationAligned(Space space) {
        super(space);
        iterator = new AtomIteratorArrayListSimple();
    }

    /**
     * Sets all spins to be aligned in the +x direction
     */
    public void initializePositions(AtomArrayList[] atomList) {
        for(int i=0; i<atomList.length; i++) {
            iterator.setList(atomList[i]);
            iterator.reset();
            while(iterator.hasNext()) {
                Vector spin = iterator.nextAtom().coord.position();
                spin.E(0.0);
                spin.setX(0,1.0);
            }
        }
    }

    private final AtomIteratorArrayListSimple iterator;
}
