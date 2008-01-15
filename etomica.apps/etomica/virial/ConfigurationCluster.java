package etomica.virial;

import etomica.atom.AtomSet;
import etomica.atom.AtomTypeMolecule;
import etomica.atom.IMolecule;
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
public class ConfigurationCluster implements Configuration, java.io.Serializable {

	/**
	 * @see etomica.config.Configuration#initializePositions(etomica.AtomIterator)
	 */
    //XXX this can't actually handle multi-atom molecules
	public void initializeCoordinates(Box box) {
        IVector dimVector = box.getSpace().makeVector();
        dimVector.E(box.getBoundary().getDimensions());
		AtomSet moleculeList = box.getMoleculeList();
		for (int i=0; i<moleculeList.getAtomCount(); i++) {
            // initialize coordinates of child atoms
		    IMolecule a = (IMolecule)moleculeList.getAtom(i);
            Conformation config = ((AtomTypeMolecule)a.getType()).getConformation();
            config.initializePositions(a.getChildList());
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
