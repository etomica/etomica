/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations.theta;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;
import etomica.virial.BoxCluster;

/**
 * Monte Carlo move for Mayer sampling that explore bond angle degrees of
 * freedom.  One bond angle is varied in each molecule.  Atoms on each side of
 * the bond are rotated around the middle atom of the bond angle with the
 * overall molecule COM fixed.
 *
 * Originally developed by Arpit for alkanes.
 */
public class MCMoveClusterShuffle extends MCMoveBoxStep {
    private final PotentialCompute potential;
    protected final IRandom random;
    protected final Space space;
    protected int numMoved;
    protected int start;
    protected Vector[] oldPosition;
    protected int[] seq;
    protected int iMolecule;
    double uOld = 0;
    double uNew = 0;
    double wOld = 0;
    double wNew = 0;

    public MCMoveClusterShuffle(PotentialCompute potentialCompute, Space space, IRandom random) {
        super();
        this.potential = potentialCompute;
        this.space = space;
        this.random = random;
        setStepSizeMin(2);
    }

    public void setBox(Box p) {
        super.setBox(p);
        int n = p.getMoleculeList().get(0).getChildList().size();
        oldPosition = space.makeVectorArray(n);
        seq = new int[n-2];
        setStepSizeMax(n-2);
    }

    @Override
    public double energyChange() {
        return 0;
    }

    @Override
    public boolean doTrial() {
        numMoved = ((int)stepSize) + (random.nextDouble() < stepSize-(int)stepSize ? 1 : 0);
        uOld = potential.computeAll(false);
        wNew = wOld = 1;
        if (box instanceof BoxCluster) wOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        IMoleculeList moleculeList = box.getMoleculeList();
        iMolecule = random.nextInt(moleculeList.size());

        IMolecule molecule = moleculeList.get(iMolecule);
        IAtomList atoms = molecule.getChildList();
        // setup
        numMoved = Math.min(numMoved, atoms.size()-1);
        start = 1 + random.nextInt(atoms.size()-numMoved-1);
        for (int j=0; j<numMoved; j++) {
            seq[j] = start+j;
            oldPosition[j].Ev1Mv2(atoms.get(start+j).getPosition(), atoms.get(start+j-1).getPosition());
        }

        Vector shift = space.makeVector();

        // actually shuffle sequence
        for (int j=0; j<numMoved; j++) {
            int s = j + random.nextInt(numMoved-j);
            int ss = seq[s];
            seq[s] = seq[j];
            seq[j] = ss;
            shift.PE(atoms.get(start+j).getPosition());
            atoms.get(start+j).getPosition().Ev1Pv2(atoms.get(start+j-1).getPosition(), oldPosition[ss-start]);
            shift.ME(atoms.get(start+j).getPosition());
        }

        if (iMolecule > 0 || true) {
            shift.TE(1.0/atoms.size());
            for (IAtom a : atoms) {
                a.getPosition().PE(shift);
            }
        }

        uNew = potential.computeAll(false);
        if (box instanceof BoxCluster) {
            ((BoxCluster)box).trialNotify();
            wNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        }
        return true;
    }

    @Override
    public double getChi(double temperature) {
        return (wOld == 0 ? 1 : wNew / wOld) * Math.exp(-(uNew - uOld) / temperature);
    }

    @Override
    public void acceptNotify() {
        if (box instanceof BoxCluster) ((BoxCluster)box).acceptNotify();
    }

    @Override
    public void rejectNotify() {
        IMoleculeList moleculeList = box.getMoleculeList();
        IAtomList atoms = moleculeList.get(iMolecule).getChildList();
        Vector shift = space.makeVector();
        for (int j=0; j<numMoved; j++) {
            shift.PE(atoms.get(start+j).getPosition());
            atoms.get(start+j).getPosition().Ev1Pv2(atoms.get(start+j-1).getPosition(), oldPosition[j]);
            shift.ME(atoms.get(start+j).getPosition());
        }

        if (iMolecule > 0 || true) {
            shift.TE(1.0/atoms.size());
            for (IAtom a : atoms) {
                a.getPosition().PE(shift);
            }
        }

        if (box instanceof BoxCluster) ((BoxCluster)box).rejectNotify();
    }
}
