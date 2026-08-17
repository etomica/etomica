/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.util.random.IRandom;
import etomica.space.Space;

public class ConfigurationClusterMove extends ConfigurationCluster {
	protected final IRandom random;
	protected final double distance;
	protected final MCMoveBox[] moves;
    public ConfigurationClusterMove(Space _space, IRandom random) {
        this(_space, random, 2);
    }
    
	public ConfigurationClusterMove(Space _space, IRandom random, double distance) {
		this(_space, random, distance, new MCMoveBox[0]);
	}
	public ConfigurationClusterMove(Space _space, IRandom random, double distance, MCMoveBox[] moves){
		super(_space);
		this.random = random;
		this.distance = distance;
		this.moves = moves;
	}

	public void initializeCoordinates(Box box) {
		super.initializeCoordinates(box);
		BoxCluster clusterBox =(BoxCluster) box;
		while (clusterBox.getSampleCluster().value(clusterBox) == 0) {
    		IAtomList list = box.getLeafList();
    		for (int i = 1; i<list.size(); i++){
    			list.get(i).getPosition().setRandomInSphere(random);
    			list.get(i).getPosition().TE(distance);
    		}
            clusterBox.trialNotify();
            clusterBox.acceptNotify();
			if(clusterBox.getSampleCluster().value(clusterBox) == 0 && moves.length > 0){
				MCMoveBox m = moves[random.nextInt(moves.length)];
				if(m.doTrial()){
					m.acceptNotify();
				}
			}
        }
	}
}
