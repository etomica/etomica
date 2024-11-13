/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.molecule.CenterOfMass;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

public class ConfigurationClusterMoveMolecule extends ConfigurationCluster {
	protected final IRandom random;
	protected final double distance;
	protected final MCMoveBox[] moves;
    public ConfigurationClusterMoveMolecule(Space _space, IRandom random) {
        this(_space, random, 2);
    }

	public ConfigurationClusterMoveMolecule(Space _space, IRandom random, double distance) {
		this(_space, random, distance, new MCMoveBox[0]);
	}
	public ConfigurationClusterMoveMolecule(Space _space, IRandom random, double distance, MCMoveBox[] moves){
		super(_space);
		this.random = random;
		this.distance = distance;
		this.moves = moves;
	}

	public void initializeCoordinates(Box box) {
		super.initializeCoordinates(box);
		BoxCluster clusterBox =(BoxCluster) box;
		Vector translationVector = clusterBox.getSpace().makeVector();
		while (clusterBox.getSampleCluster().value(clusterBox) == 0) {
    		IMoleculeList list = box.getMoleculeList();
    		for (int i = 1; i<list.size() - 1; i++){
				IMolecule molecule = list.get(i);
				Vector com = CenterOfMass.position(clusterBox, molecule);
				translationVector.setRandomInSphere(random);
				translationVector.TE(distance);
				translationVector.ME(com);
				molecule.getChildList().forEach(atom -> {

					atom.getPosition().PE(translationVector);
				});
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
