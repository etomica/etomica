/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.mcmove;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.util.collections.IntArrayList;
import etomica.util.random.IRandom;
import etomica.virial.BoxCluster;

/**
 * Monte Carlo move for Mayer sampling that explore bond length degrees of
 * freedom.  One bond length is varied in each molecule.  Atoms on each side of
 * the bond are translated with the overall molecule COM fixed.
 *
 * Originally developed by Arpit for alkanes.
 */
public class MCMoveClusterStretch extends MCMoveBoxStep {
    private final PotentialCompute potential;
    protected final IRandom random;
    protected final Space space;
    protected ISpecies species;
    protected double step = 0;
    protected Vector[][] position = null;
    protected final IntArrayList [] bonding;
    double uOld = 0;
    double uNew = 0;
    double wOld = 0;
    double wNew = 0;
    int [] modified;
    int modifiedIndex = 0;
    int b = 0;

    public MCMoveClusterStretch(PotentialCompute potentialCompute, Space space, IntArrayList[] bonding, IRandom random, double stepSize) {
        super();
        this.potential = potentialCompute;
        this.space = space;
        this.random = random;
        this.bonding = bonding;
        this.stepSize = stepSize;
        modified = new int[bonding.length];
    }

    public void setBox(Box p) {
        super.setBox(p);
        position = new Vector[p.getMoleculeList().size()][p.getMoleculeList().get(0).getChildList().size()];
    }

    @Override
    public double energyChange() {
        return 0;
    }

    @Override
    public boolean doTrial() {
        uOld = potential.computeAll(false);
        wOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        IMoleculeList moleculeList = box.getMoleculeList();
        for(int i = 0; i<moleculeList.size(); i++) {
            if (species != null && moleculeList.get(i).getType() != species) {
                continue;
            }
            IMolecule molecule = moleculeList.get(i);
            IAtomList atoms = molecule.getChildList();
            for(int j = 0; j < molecule.getChildList().size(); j++) {
                position[i][j] = space.makeVector();
                position[i][j].E(molecule.getChildList().get(j).getPosition());
            }
            modifiedIndex = 0;
            step = stepSize * (random.nextDouble() - 0.5);
            do{
                b = random.nextInt(bonding.length);
            }while (bonding[b].size() < 1);
            modified[modifiedIndex]=b;
            ++modifiedIndex;
            int a = random.nextInt(bonding[b].size());
            a = bonding[b].getInt(a);
            modified[modifiedIndex] = a;
            ++modifiedIndex;
            Vector dr = space.makeVector();
            dr.Ev1Mv2(atoms.get(b).getPosition(), atoms.get(a).getPosition());
            dr.TE(step / Math.sqrt(dr.squared()));
            Vector shift = space.makeVector();
            double m = 0;
            m += transform(dr, atoms.get(b), shift);
            m += transformBondedAtoms(dr, b, atoms, shift);
            dr.TE(-1);
            m += transform(dr, atoms.get(a), shift);
            m += transformBondedAtoms(dr, a, atoms, shift);
            shift.TE(-1.0/m);
            for (IAtom aa : atoms) {
                aa.getPosition().PE(shift);
            }
        }
        ((BoxCluster)box).trialNotify();
        uNew = potential.computeAll(false);
        wNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        return true;
    }

    protected double transformBondedAtoms(Vector dr, int index, IAtomList atoms, Vector shift){
        double m = 0;
        for(int k = 0; k < bonding[index].size(); k++){
            boolean rotated = false;
            for (int l = 0; l < modifiedIndex; l++) {
                if (bonding[index].getInt(k) == modified[l]) {
                    rotated = true;
                    break;
                }
            }
            if (!rotated) {
                m += transform(dr, atoms.get(bonding[index].getInt(k)), shift);
                modified[modifiedIndex] = bonding[index].getInt(k);
                ++modifiedIndex;
                transformBondedAtoms(dr, bonding[index].getInt(k), atoms, shift);
            }
        }
        return m;
    }

    protected double transform(Vector dr, IAtom a, Vector shift) {
        a.getPosition().PE(dr);
        double m = a.getType().getMass();
        shift.PEa1Tv1(m, dr);
        return m;
    }

    @Override
    public double getChi(double temperature) {
        return (wOld == 0 ? 1 : wNew / wOld) * Math.exp(-(uNew - uOld) / temperature);
    }

    @Override
    public void acceptNotify() {
        ((BoxCluster)box).acceptNotify();
    }

    @Override
    public void rejectNotify() {
        IMoleculeList moleculeList = box.getMoleculeList();
        for(int i = 0; i<box.getMoleculeList().size(); i++) {
            for(int j = 0; j < moleculeList.get(i).getChildList().size(); j++) {
                moleculeList.get(i).getChildList().get(j).getPosition().E(position[i][j]);
            }
        }
        ((BoxCluster)box).rejectNotify();
    }
}
