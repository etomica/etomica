/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculePositionCOM;
import etomica.space.Space;
import etomica.util.random.IRandom;

/**
 * Class that samples molecule positions based on a chain of hard spheres of
 * diameter sigma.  The sequence of spheres is chosen randomly.
 *
 * @author Andrew
 */
public class MCMoveClusterMoleculeHSChain extends MCMoveMolecule {

    protected final double sigma;
    protected int[] seq;
    protected MoleculePositionCOM moleculeCOM ;

    public MCMoveClusterMoleculeHSChain(IRandom random, Space _space, double sigma) {
        super(null,random, _space,1, 1);
        this.sigma = sigma;
        moleculeCOM = new MoleculePositionCOM(space);
    }

    public boolean doTrial() {

        IMoleculeList moleculeList = box.getMoleculeList();
        int n = moleculeList.size();
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
        IMolecule molecule = moleculeList.get(seq[0]);
        groupTranslationVector.E(moleculeCOM.position(molecule));

        groupTranslationVector.TE(-1);

        moveMoleculeAction.actionPerformed(molecule);

        for (int i = 1; i<n; i++) {
            groupTranslationVector.setRandomInSphere(random);

            double sig = getSigma(seq[i - 1], seq[i]);
            groupTranslationVector.TE(sig);

            IMolecule moleculePrevious = moleculeList.get(seq[i-1]);
            groupTranslationVector.PE(moleculeCOM.position(moleculePrevious));

            molecule = moleculeList.get(seq[i]);
            groupTranslationVector.ME(moleculeCOM.position(molecule));

            moveMoleculeAction.actionPerformed(molecule);
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
}
