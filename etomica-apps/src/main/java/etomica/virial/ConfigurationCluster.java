/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.box.Box;
import etomica.config.Configuration;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * @author kofke
 *
 * Generates a configuration such that the value of the box's
 * sampling cluster is positive at beta = 1.
 */
public class ConfigurationCluster implements Configuration, java.io.Serializable {

	public ConfigurationCluster(Space _space) {
		this.space = _space;
	}

	/**
	 * @see etomica.config.Configuration#initializeCoordinates
	 */
	public void initializeCoordinates(Box box) {
        Vector dimVector = space.makeVector();
        dimVector.E(box.getBoundary().getBoxSize());
		IMoleculeList moleculeList = box.getMoleculeList();
		for (int i = 0; i<moleculeList.size(); i++) {
            // initialize coordinates of child atoms
		    IMolecule a = moleculeList.get(i);
            a.getType().initializeConformation(a);
        }

        BoxCluster boxCluster = (BoxCluster)box;
        boxCluster.trialNotify();
        boxCluster.acceptNotify();
        
        // All the molecules are now at the origin.  If this isn't enough,
        // you'll need to do something more, perhaps a subclass.  Currently,
        // nothing needs that (and alkanes are unhappy with it).
	}

	protected final Space space;
    private static final long serialVersionUID = 3L;
}
