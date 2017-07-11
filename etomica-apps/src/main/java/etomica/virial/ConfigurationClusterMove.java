/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.util.random.IRandom;
import etomica.space.Space;

public class ConfigurationClusterMove extends ConfigurationCluster {

    public ConfigurationClusterMove(Space _space, IRandom random) {
        this(_space, random, 2);
    }
    
	public ConfigurationClusterMove(Space _space, IRandom random, double distance) {
		super(_space);
		this.random = random;
		this.distance = distance;
	}

	public void initializeCoordinates(Box box) {
		super.initializeCoordinates(box);
		BoxCluster clusterBox =(BoxCluster) box;
		while (clusterBox.getSampleCluster().value(clusterBox) == 0) {
    		IAtomList list = box.getLeafList();
    		for (int i=1;i<list.getAtomCount();i++){
    			list.getAtom(i).getPosition().setRandomInSphere(random);
    			list.getAtom(i).getPosition().TE(distance);
    		}
            clusterBox.trialNotify();
            clusterBox.acceptNotify();
        }
	}
   protected final IRandom random;
   protected final double distance;
}
