/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations.theta;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.molecule.CenterOfMass;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.collections.IntArrayList;
import etomica.util.random.IRandom;
import etomica.virial.BoxCluster;

/**
 * Monte Carlo move for Mayer (or Boltzmann) sampling that shuffles the bonds
 * among some group of atoms in a linear section of a polymer molecule.  The
 * endpoints of the section remain fixed.
 *
 * For nonlinear topologies, max step size needs to be specified.
 */
public class MCMoveClusterShuffle extends MCMoveBoxStep {
    private final PotentialCompute potential;
    protected final IRandom random;
    protected final Space space;
    protected int numMoved;
    protected Vector[] bondVector;
    protected int[] seq, imposedBonds;
    protected int iMolecule;
    double uOld = 0;
    double uNew = 0;
    double wOld = 0;
    double wNew = 0;
    protected boolean doLattice;
    protected IntArrayList[] bonding;

    public MCMoveClusterShuffle(PotentialCompute potentialCompute, Space space, IRandom random) {
        super();
        this.potential = potentialCompute;
        this.space = space;
        this.random = random;
        setStepSizeMin(2);
    }

    public void setBonding(IntArrayList[] bonding) {
        this.bonding = bonding;
    }

    public void setDoLattice(boolean doLattice) {
        this.doLattice = doLattice;
    }

    public void setBox(Box p) {
        super.setBox(p);
        int n = p.getMoleculeList().get(0).getChildList().size();
        bondVector = space.makeVectorArray(n);
        seq = new int[n-1];
        imposedBonds = new int[n-1];
        if (bonding == null) {
            bonding = new IntArrayList[n];
            bonding[0] = new IntArrayList(new int[]{1});
            for (int i=1; i<n-1; i++) {
                bonding[i] = new IntArrayList(new int[]{i-1,i+1});
            }
            bonding[n-1] = new IntArrayList(new int[]{n-2});
            setStepSizeMax(n-2);
        }
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
        boolean forward = false;
        int start = 1 + random.nextInt(atoms.size() - numMoved - 1);
        int actualMoved = 1;
        while (actualMoved < numMoved) {
            actualMoved = 1;
            while (bonding[start].size() != 2) {
                start = 1 + random.nextInt(atoms.size() - numMoved - 1);
            }
            seq[1] = start;
            forward = random.nextInt(2) == 0;
            seq[0] = bonding[start].getInt(forward ? 0 : 1);
            // bondVector[i] is vector forward from i
            bondVector[0].Ev1Mv2(atoms.get(seq[1]).getPosition(), atoms.get(seq[0]).getPosition());
            for (int i = 2; i <= numMoved && bonding[seq[i - 1]].size() == 2; i++) {
                seq[i] = bonding[seq[i - 1]].getInt(forward ? 1 : 0);
                bondVector[i - 1].Ev1Mv2(atoms.get(seq[i]).getPosition(), atoms.get(seq[i - 1]).getPosition());
                actualMoved++;
            }
        }

        for (int i=0; i<numMoved; i++) {
            imposedBonds[i] = i;
        }
        Vector shift = space.makeVector();

        // actually shuffle sequence
        for (int j=0; j<numMoved-1; j++) {
            int s = j + random.nextInt(numMoved-j);
            // we pick a new bond vector to impose on seq j
            int ss = imposedBonds[s];
            imposedBonds[s] = imposedBonds[j];
            imposedBonds[j] = ss;
            shift.PE(atoms.get(seq[j]).getPosition());
            atoms.get(seq[j+1]).getPosition().Ev1Pv2(atoms.get(seq[j]).getPosition(), bondVector[ss]);
            shift.ME(atoms.get(seq[j]).getPosition());
        }

        if (doLattice && iMolecule==0) {
            shift.E(CenterOfMass.position(box, molecule));
            shift.TE(-1);
            for (int k=0; k<shift.getD(); k++) {
                shift.setX(k, Math.floor(shift.getX(k)));
            }
        }
        else if (!doLattice) {
            shift.TE(1.0/atoms.size());
        }
        for (IAtom a : atoms) {
            a.getPosition().PE(shift);
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
            shift.PE(atoms.get(seq[j]).getPosition());
            atoms.get(seq[j+1]).getPosition().Ev1Pv2(atoms.get(seq[j]).getPosition(), bondVector[j]);
            shift.ME(atoms.get(seq[j]).getPosition());
        }

        if (doLattice && iMolecule==0) {
            shift.E(CenterOfMass.position(box, moleculeList.get(iMolecule)));
            shift.TE(-1);
            for (int k=0; k<shift.getD(); k++) {
                shift.setX(k, Math.floor(shift.getX(k)));
            }
        }
        else if (!doLattice) {
            shift.TE(1.0/atoms.size());
        }
        for (IAtom a : atoms) {
            a.getPosition().PE(shift);
        }

        if (box instanceof BoxCluster) ((BoxCluster)box).rejectNotify();
    }
}
