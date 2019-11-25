/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.IAtomList;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.space.Space;
import etomica.util.random.IRandom;

/**
 * Class that samples molecule positions by randomly placing hard spheres
 * in the box.
 *
 * @author Arpit Bansal
 */
public class MCMoveClusterAtomInBox extends MCMoveAtom {

    public MCMoveClusterAtomInBox(IRandom random, Space _space){//}, double sigma) {
        super(random, null, _space);
        //this.sigma = sigma;
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
}
