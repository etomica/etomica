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
import etomica.space3d.RotationTensor3D;
import etomica.species.ISpecies;
import etomica.util.collections.IntArrayList;
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
public class MCMoveClusterAngleGeneral extends MCMoveBoxStep {
    private final PotentialCompute potential;
    protected final IRandom random;
    protected final Space space;
    protected final ISpecies species;
    protected double dt = 0;
    protected Vector[] position = null;
    protected final int[][] triplets;
    protected int iMolecule;
    protected final IntArrayList[] bonding;
    double uOld = 0;
    double uNew = 0;
    double wOld = 0;
    double wNew = 0;
    final int [] modified;
    int modifiedIndex = 0;
    int b = 0;
    protected int[] constraintMap;

    public MCMoveClusterAngleGeneral(PotentialCompute potentialCompute, Space space, ISpecies species, IntArrayList[] bonding, int[][] triplets, IRandom random, double stepSize) {
        super();
        this.potential = potentialCompute;
        this.space = space;
        this.random = random;
        this.triplets = triplets;
        this.stepSize = stepSize;
        this.species = species;
        this.bonding = bonding;
        modified = new int[species.getLeafAtomCount()];
        setStepSizeMax(Math.PI/2);
    }


    public void setBox(Box p) {
        super.setBox(p);
        position = space.makeVectorArray(p.getMoleculeList(species).get(0).getChildList().size());
    }


    public void setConstraintMap(int[] newConstraintMap) {
        constraintMap = newConstraintMap;
    }

    @Override
    public double energyChange() {
        return 0;
    }

    @Override
    public boolean doTrial() {
        uOld = potential.computeAll(false);
        wNew = wOld = 1;
        if (box instanceof BoxCluster) wOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        IMoleculeList moleculeList = box.getMoleculeList();
        iMolecule = random.nextInt(moleculeList.size());
        while (species != null && moleculeList.get(iMolecule).getType() != species) {
            iMolecule = random.nextInt(moleculeList.size());
        }

        IMolecule molecule = moleculeList.get(iMolecule);
        IAtomList atoms = molecule.getChildList();
        for(int j = 0; j < molecule.getChildList().size(); j++) {
            position[j].E(atoms.get(j).getPosition());
        }
        modifiedIndex = 0;

        int d = random.nextInt(triplets.length);
        b = triplets[d][1];
        modified[modifiedIndex]=b;
        ++modifiedIndex;
        Vector axis = space.makeVector();
        Vector temp = space.makeVector();
        axis.Ev1Mv2(atoms.get(triplets[d][0]).getPosition(), atoms.get(b).getPosition());
        temp.Ev1Mv2(atoms.get(triplets[d][2]).getPosition(), atoms.get(b).getPosition());
        axis.XE(temp);
        axis.normalize();

        dt = 2 * stepSize * (random.nextDouble() - 0.5);

        RotationTensor3D rotationTensor = new RotationTensor3D();
        rotationTensor.setRotationAxis(axis, dt);
        Vector shift = space.makeVector();
        transformBondedAtoms(rotationTensor, triplets[d][0], atoms, shift);
        rotationTensor.setRotationAxis(axis, -dt);

        transformBondedAtoms(rotationTensor, triplets[d][2], atoms, shift);

        if (iMolecule==0 || (constraintMap != null && constraintMap[iMolecule] == 0)) {
            double mt = 0;
            for (IAtom aa : atoms) {
                mt += aa.getType().getMass();
            }
            shift.TE(-1.0 / mt);
            for (IAtom aa : atoms) {
                aa.getPosition().PE(shift);
            }
        }

        if (box instanceof BoxCluster) {
            ((BoxCluster)box).trialNotify();
            wNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        }
        uNew = potential.computeAll(false);

        return true;
    }

    protected void transformBondedAtoms(RotationTensor3D rotationTensor3D, int index, IAtomList atoms, Vector shift){
        for(int k = 0; k < bonding[index].size(); k++){
            boolean rotated = false;
            for (int l = 0; l < modifiedIndex; l++) {
                if (bonding[index].getInt(k) == modified[l]) {
                    rotated = true;
                    break;
                }
            }
            if (!rotated) {
                transform(rotationTensor3D, bonding[index].getInt(k), atoms, shift);
                modified[modifiedIndex] = bonding[index].getInt(k);
                ++modifiedIndex;
                transformBondedAtoms(rotationTensor3D, bonding[index].getInt(k), atoms, shift);
            }
        }
    }

    protected void transform(RotationTensor3D rotationTensor3D, int index, IAtomList atoms, Vector shift) {
        Vector r = space.makeVector();
        IAtom a = atoms.get(index);
        double m = a.getType().getMass();
        Vector p = a.getPosition();
        shift.PEa1Tv1(-m, p);
        r.Ev1Mv2(p, atoms.get(b).getPosition());
        rotationTensor3D.transform(r);
        r.PE(atoms.get(b).getPosition());
        p.E(r);
        shift.PEa1Tv1(+m, p);
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
        for(int j = 0; j < moleculeList.get(iMolecule).getChildList().size(); j++) {
            moleculeList.get(iMolecule).getChildList().get(j).getPosition().E(position[j]);
        }
        if (box instanceof BoxCluster) ((BoxCluster)box).rejectNotify();
    }
}
