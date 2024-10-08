/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.mcmove;

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
public class MCMoveClusterAngle extends MCMoveBoxStep {
    private final PotentialCompute potential;
    protected final IRandom random;
    protected final Space space;
    protected ISpecies species;
    protected double dt = 0;
    protected Vector[] position = null;
    protected final IntArrayList [] bonding;
    protected int iMolecule;
    double uOld = 0;
    double uNew = 0;
    double wOld = 0;
    double wNew = 0;
    int [] modified;
    int modifiedIndex = 0;
    int b = 0;
    protected int start, stop;
    protected boolean doLattice;
    protected int[] constraintMap;

    public MCMoveClusterAngle(PotentialCompute potentialCompute, Space space, IntArrayList[] bonding, IRandom random, double stepSize) {
        super();
        this.potential = potentialCompute;
        this.space = space;
        this.random = random;
        this.bonding = bonding;
        this.stepSize = stepSize;
        modified = new int[bonding.length];
        setStepSizeMax(Math.PI/2);
        start = 0;
        stop = Integer.MAX_VALUE;
    }

    public void setDoLattice(boolean doLattice) {
        this.doLattice = doLattice;
    }

    public void setBox(Box p) {
        super.setBox(p);
        position = space.makeVectorArray(p.getMoleculeList().get(0).getChildList().size());
    }

    public void setAtomRange(int start, int stop) {
        this.start = start;
        this.stop = stop;
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
        if (molecule.getChildList().size() < 3) return false;
        modifiedIndex = 0;
        int d = 0;
        do{
            b = random.nextInt(bonding.length);
            d = Math.min(b, bonding.length-1-b);
        }while (bonding[b].size() < 2 || (d < start || d >= stop));
        modified[modifiedIndex]=b;
        ++modifiedIndex;
        int a = random.nextInt(bonding[b].size());
        a = bonding[b].getInt(a);
        modified[modifiedIndex] = a;
        ++modifiedIndex;
        Vector r = space.makeVector();
        r.Ev1Mv2(atoms.get(b).getPosition(), atoms.get(a).getPosition());
        Vector axis = space.makeVector();
        if (doLattice) {
            dt = Math.PI/2;
            do {
                axis.E(0);
                axis.setX(random.nextInt(axis.getD()), random.nextInt(2) * 2 - 1);
            }
            while (Math.abs(axis.dot(r))/Math.sqrt(r.squared()) > 0.9);
        }
        else {
            dt = 2 * stepSize * (random.nextDouble() - 0.5);
            axis.setRandomSphere(random);
            Vector projection = space.makeVector();
            projection.Ea1Tv1(axis.dot(r)/r.squared(), r);
            axis.ME(projection);
            axis.normalize();
        }
        RotationTensor3D rotationTensor = new RotationTensor3D();
        rotationTensor.setRotationAxis(axis, dt);
        Vector shift = space.makeVector();
        transform(rotationTensor, a, atoms, shift);
        transformBondedAtoms(rotationTensor, a, atoms, shift);

        if (iMolecule==0 || (constraintMap != null && constraintMap[iMolecule] == 0)) {
            if (doLattice) {
                shift.E(CenterOfMass.position(box, molecule));
                shift.TE(-1);
                for (int i=0; i<shift.getD(); i++) {
                    shift.setX(i, Math.floor(shift.getX(i)));
                }
            }
            else {
                double mt = 0;
                for (IAtom aa : atoms) {
                    mt += aa.getType().getMass();
                }
                shift.TE(-1.0 / mt);
            }
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
