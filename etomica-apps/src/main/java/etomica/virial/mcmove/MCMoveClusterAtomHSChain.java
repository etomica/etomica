/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.mcmove;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.space.Vector;
import etomica.util.random.IRandom;
import etomica.virial.BoxCluster;

/**
 * Class that samples molecule positions based on a chain of hard spheres of
 * diameter sigma.  The sequence of spheres is chosen randomly.
 *
 * @author Andrew
 */
public class MCMoveClusterAtomHSChain extends MCMoveBoxStep {

    protected final IRandom random;
    protected final double sigma;
    protected int[] seq;
    protected boolean forceInBox;

    public MCMoveClusterAtomHSChain(IRandom random, Box box, double sigma) {
        super();
        this.random = random;
        this.sigma = sigma;
        setBox(box);
    }

    public void setForceInBox(boolean forceInBox) {
        this.forceInBox = forceInBox;
    }

    public boolean doTrial() {

        IAtomList leafAtoms = box.getLeafList();
        int n = leafAtoms.size();
        if (seq == null) {
            seq = new int[n];
        }
        for (int i = 0; i<n; i++) {
            seq[i] = i;
        }
        for (int i = 0; i<n; i++) {
            int j = i+random.nextInt(n-i);
            int k = seq[j];
            seq[j] = seq[i];
            seq[i] = k;
        }
        leafAtoms.get(seq[0]).getPosition().E(0);

        for (int i = 1; i<n; i++) {
            Vector pos = leafAtoms.get(seq[i]).getPosition();
            boolean insideBox = false;
            while(!insideBox) {
                pos.setRandomInSphere(random);
                double sig = getSigma(seq[i - 1], seq[i]);
                if (sig < 0) {
                    // we want to force the position to be in the well (between 1 and sigma)
                    sig = -sig;
                    while (pos.squared() < 1 / (sig * sig)) {
                        pos.setRandomInSphere(random);
                    }
                }
                pos.TE(sig);
                insideBox = true;
                if (forceInBox) {
                    for (int j = 0; j < pos.getD(); j++) {
                        if (Math.abs(pos.getX(j)) > box.getBoundary().getBoxSize().getX(j) / 2) {
                            insideBox = false;
                            break;
                        }
                    }
                }
                if (insideBox) {
                    pos.PE(leafAtoms.get(seq[i - 1]).getPosition());
                }
            }
        }

        ((BoxCluster)box).trialNotify();
        return true;
    }

    protected double getSigma(int i, int j) {
        return sigma;
    }

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
