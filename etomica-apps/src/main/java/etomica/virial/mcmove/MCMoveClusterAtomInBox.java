/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.mcmove;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.util.random.IRandom;
import etomica.virial.BoxCluster;

/**
 * Class that samples molecule positions by randomly placing hard spheres
 * in the box.
 *
 * @author Arpit Bansal
 */
public class MCMoveClusterAtomInBox extends MCMoveBox {

    protected final IRandom random;

    public MCMoveClusterAtomInBox(IRandom random, Box box){
        super();
        this.random = random;
        setBox(box);
    }

    public boolean doTrial() {

        IAtomList leafAtoms = box.getLeafList();
        int n = box.getLeafList().size();
        for (int i=0; i<n; i++){
            leafAtoms.get(i).getPosition().setRandomCube(random);
            leafAtoms.get(i).getPosition().TE(box.getBoundary().getBoxSize());
        }
        ((BoxCluster)box).trialNotify();
        return true;
    }

    // override this to do a mixture

    public double getChi(double temperature) {
        return 1;
    }

    public void rejectNotify() {
        throw new RuntimeException("nope");
    }

    public void acceptNotify() {
        ((BoxCluster)box).acceptNotify();
    }

    @Override
    public double energyChange() {
        return 0;
    }
}
