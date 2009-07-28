package etomica.virial;

import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.IVectorMutable;
import etomica.config.Configuration;
import etomica.config.IConformation;
import etomica.space.ISpace;

/**
 * @author kofke
 *
 * Generates a configuration such that the value of the box's
 * sampling cluster is positive at beta = 1.
 */
public class ConfigurationCluster implements Configuration, java.io.Serializable {

	public ConfigurationCluster(ISpace _space) {
		this.space = _space;
	}

	/**
	 * @see etomica.config.Configuration#initializePositions(etomica.AtomIterator)
	 */
    //XXX this can't actually handle multi-atom molecules
	public void initializeCoordinates(IBox box) {
        IVectorMutable dimVector = space.makeVector();
        dimVector.E(box.getBoundary().getBoxSize());
		IMoleculeList moleculeList = box.getMoleculeList();
		for (int i=0; i<moleculeList.getMoleculeCount(); i++) {
            // initialize coordinates of child atoms
		    IMolecule a = moleculeList.getMolecule(i);
            a.getType().initializeConformation(a);
        }

        BoxCluster boxCluster = (BoxCluster)box;
        boxCluster.trialNotify();
        boxCluster.acceptNotify();
        
        // All the molecules are now at the origin.  If this isn't enough,
        // you'll need to do something more, perhaps a subclass.  Currently,
        // nothing needs that (and alkanes are unhappy with it).
	}

	private final ISpace space;
    private static final long serialVersionUID = 3L;
}
