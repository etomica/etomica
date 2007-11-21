package etomica.virial;

import etomica.atom.AtomTypeGroup;
import etomica.atom.IAtom;
import etomica.atom.IAtomGroup;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.config.Conformation;
import etomica.space.IVector;

/**
 * @author kofke
 *
 * Generates a configuration such that the value of the box's
 * sampling cluster is positive at beta = 1.
 */
public class ConfigurationCluster extends Configuration {

	/**
	 * @see etomica.config.Configuration#initializePositions(etomica.AtomIterator)
	 */
	public void initializeCoordinates(Box box) {
        IVector dimVector = box.getSpace().makeVector();
        dimVector.E(box.getBoundary().getDimensions());
        AtomIteratorAllMolecules iterator = new AtomIteratorAllMolecules(box);
		iterator.reset();
        for (IAtom a = iterator.nextAtom(); a != null;
             a = iterator.nextAtom()) {
            if (a instanceof IAtomGroup) {
                // initialize coordinates of child atoms
                Conformation config = ((AtomTypeGroup)a.getType()).getConformation();
                config.initializePositions(((IAtomGroup)a).getChildList());
            }
        }
        BoxCluster boxCluster = (BoxCluster)box;
        boxCluster.trialNotify();
        boxCluster.acceptNotify();
        
        // All the molecules are now at the origin.  If this isn't enough,
        // you'll need to do something more, perhaps a subclass.  Currently,
        // nothing needs that (and alkanes are unhappy with it).
	}

    private static final long serialVersionUID = 3L;
}
