/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.mcmove;

import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.molecule.CenterOfMass;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.space.Vector;
import etomica.util.random.IRandom;
import etomica.virial.BoxCluster;

/**
 * Class that samples molecule positions based on a chain of hard spheres of
 * diameter sigma.  The sequence of spheres is chosen randomly.
 *
 * @author Andrew
 */
public class MCMoveClusterMoleculeHSChain extends MCMoveBox {

    protected final IRandom random;
    protected final double sigma;
    protected int[] seq;
    protected int[] constraintMap;
    protected final Vector translationVector;

    public MCMoveClusterMoleculeHSChain(IRandom random, Box box, double sigma) {
        super();
        this.random = random;
        this.sigma = sigma;
        translationVector = box.getSpace().makeVector();
        setBox(box);
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

    @Override
    public double energyChange() {
        return 0;
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
        translationVector.E(CenterOfMass.position(box, molecule));

        translationVector.TE(-1);

        molecule.getChildList().forEach(atom -> {
            atom.getPosition().PE(translationVector);
        });
        for (int j=0; j<n; j++) {
            if (constraintMap[j] == prev && j != prev) {
                translationVector.E(0);
                translateFrom(prev, j);
            }
        }

        for (int i = start; i<n; i++) {
            if (constraintMap[seq[i]] != seq[i]) continue;
            translationVector.setRandomInSphere(random);
//            translationVector.setX(1, 0);
//            translationVector.setX(2, 0);

            double sig = getSigma(prev, seq[i]);
            translationVector.TE(sig);
            translateFrom(prev, seq[i]);

            for (int j=0; j<n; j++) {
                if (constraintMap[j] == seq[i] && j != seq[i]) {
                    translationVector.E(0);
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
        translationVector.PE(CenterOfMass.position(box, moleculePrevious));

        IMolecule molecule = moleculeList.get(current);
        translationVector.ME(CenterOfMass.position(box, molecule));

        molecule.getChildList().forEach(atom -> {
            atom.getPosition().PE(translationVector);
        });
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
