/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.mcmove;

import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculePositionCOM;
import etomica.space.Space;
import etomica.util.random.IRandom;
import etomica.virial.BoxCluster;

/**
 * Class that samples molecule positions based on a chain of hard spheres of
 * diameter sigma.  The sequence of spheres is chosen randomly.
 *
 * @author Andrew
 */
public class MCMoveClusterMoleculeHSChain extends MCMoveMolecule {

    protected final double sigma;
    protected int[] seq;
    protected MoleculePositionCOM moleculeCOM;
    protected int[] constraintMap;

    public MCMoveClusterMoleculeHSChain(IRandom random, Space _space, double sigma) {
        super(null,random, _space,1, 1);
        this.sigma = sigma;
        moleculeCOM = new MoleculePositionCOM(space);
    }

    public void setBox(Box p) {
        super.setBox(p);
        if (constraintMap == null) {
            constraintMap = new int[box.getMoleculeList().size()];
            for (int i=0; i<constraintMap.length; i++) {
                constraintMap[i] = i;
            }
        }
    }


    public void setConstraintMap(int[] newConstraintMap) {
        constraintMap = newConstraintMap;
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

        int prev = seq[0];
        int start = 1;
        if (constraintMap[prev] != prev) {
            prev = seq[1];
            start = 2;
        }
        IMolecule molecule = moleculeList.get(prev);
        groupTranslationVector.E(moleculeCOM.position(molecule));

        groupTranslationVector.TE(-1);

        moveMoleculeAction.actionPerformed(molecule);

        for (int j=0; j<n; j++) {
            if (constraintMap[j] == prev && j != prev) {
                groupTranslationVector.E(0);
                translateFrom(prev, j);
            }
        }

        for (int i = start; i<n; i++) {
            if (constraintMap[seq[i]] != seq[i]) continue;
            groupTranslationVector.setRandomInSphere(random);

            double sig = getSigma(prev, seq[i]);
            groupTranslationVector.TE(sig);
            translateFrom(prev, seq[i]);

            for (int j=0; j<n; j++) {
                if (constraintMap[j] == seq[i] && j != seq[i]) {
                    groupTranslationVector.E(0);
                    translateFrom(seq[i], j);
                }
            }
            prev = seq[i];
        }

        ((BoxCluster)box).trialNotify();
        return true;
    }

    protected void translateFrom(int prev, int current) {
        IMoleculeList moleculeList = box.getMoleculeList();
        IMolecule moleculePrevious = moleculeList.get(prev);
        groupTranslationVector.PE(moleculeCOM.position(moleculePrevious));

        molecule = moleculeList.get(current);
        groupTranslationVector.ME(moleculeCOM.position(molecule));

        moveMoleculeAction.actionPerformed(molecule);
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
