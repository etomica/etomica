/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations.theta;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.molecule.CenterOfMass;
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
public class MCMoveClusterReptate extends MCMoveBox {
    private final PotentialCompute potential;
    protected final IRandom random;
    protected final Space space;
    protected boolean[] forward;
    protected Vector[] oldPosition;
    double uOld = 0;
    double uNew = 0;
    double wOld = 0;
    double wNew = 0;
    protected boolean doLattice;

    public MCMoveClusterReptate(PotentialCompute potentialCompute, Space space, IRandom random) {
        super();
        this.potential = potentialCompute;
        this.space = space;
        this.random = random;
    }

    public void setDoLattice(boolean doLattice) {
        this.doLattice = doLattice;
    }

    public void setBox(Box p) {
        super.setBox(p);
        oldPosition = space.makeVectorArray(p.getMoleculeList().size());
        forward = new boolean[p.getMoleculeList().size()];
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
        for(int i = 0; i<moleculeList.size(); i++) {
            IMolecule molecule = moleculeList.get(i);
            IAtomList atoms = molecule.getChildList();
            forward[i] = random.nextInt(2) == 0;
            Vector drNew = space.makeVector();
            if (doLattice) {
                boolean allZero = false;
                do {
                    for (int j = 0; j < drNew.getD(); j++) {
                        int dx = random.nextInt(3) - 1;
                        drNew.setX(j, dx);
                    }
                }
                while (drNew.isZero());
            }
            else {
                drNew.setRandomSphere(random);
            }
            Vector shift = space.makeVector();
            if (forward[i]) {
                shift.E(atoms.get(0).getPosition());
                oldPosition[i].Ev1Mv2(atoms.get(0).getPosition(), atoms.get(1).getPosition());
                goForward(atoms, drNew);
                shift.ME(atoms.get(atoms.size()-1).getPosition());
            }
            else {
                shift.E(atoms.get(atoms.size()-1).getPosition());
                oldPosition[i].Ev1Mv2(atoms.get(atoms.size()-1).getPosition(), atoms.get(atoms.size()-2).getPosition());
                goBackward(atoms, drNew);
                shift.ME(atoms.get(0).getPosition());
            }
            if (i==0) {
                if (doLattice) shiftLatticeCOM(molecule, shift);
                else shiftCOM(atoms, shift);
            }
        }
        uNew = potential.computeAll(false);
        if (box instanceof BoxCluster) {
            ((BoxCluster)box).trialNotify();
            wNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        }
        return true;
    }

    protected void shiftLatticeCOM(IMolecule molecule, Vector shift) {
        shift.E(CenterOfMass.position(box, molecule));
        shift.TE(-1);
        for (int k=0; k<shift.getD(); k++) {
            shift.setX(k, Math.floor(shift.getX(k)));
        }
        for (IAtom aa : molecule.getChildList()) {
            aa.getPosition().PE(shift);
        }
    }

    protected void shiftCOM(IAtomList atoms, Vector shift) {
        double mt = atoms.get(0).getType().getMass() * atoms.size();
        shift.TE(1.0/mt);
        for (IAtom aa : atoms) {
            aa.getPosition().PE(shift);
        }
    }

    protected void goForward(IAtomList atoms, Vector drNew) {
        for (int i=0; i<atoms.size()-1; i++) {
            atoms.get(i).getPosition().E(atoms.get(i+1).getPosition());
        }
        atoms.get(atoms.size()-1).getPosition().Ev1Pv2(atoms.get(atoms.size()-2).getPosition(), drNew);
    }

    protected void goBackward(IAtomList atoms, Vector drNew) {
        for (int i=atoms.size()-1; i>0; i--) {
            atoms.get(i).getPosition().E(atoms.get(i-1).getPosition());
        }
        atoms.get(0).getPosition().Ev1Pv2(atoms.get(1).getPosition(), drNew);

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
        for(int i = 0; i<box.getMoleculeList().size(); i++) {
            IAtomList atoms = moleculeList.get(i).getChildList();
            Vector rNew = space.makeVector();
            Vector shift = space.makeVector();
            if (forward[i]) {
                shift.E(atoms.get(atoms.size()-1).getPosition());
                rNew.Ev1Mv2(atoms.get(atoms.size()-1).getPosition(), atoms.get(atoms.size()-2).getPosition());
                goBackward(atoms, oldPosition[i]);
                shift.ME(atoms.get(0).getPosition());
            }
            else {
                shift.E(atoms.get(0).getPosition());
                rNew.Ev1Mv2(atoms.get(0).getPosition(), atoms.get(1).getPosition());
                goForward(atoms, oldPosition[i]);
                shift.ME(atoms.get(atoms.size()-1).getPosition());
            }
            if (i==0) {
                if (doLattice) shiftLatticeCOM(moleculeList.get(i), shift);
                else shiftCOM(atoms, shift);
            }
        }
        if (box instanceof BoxCluster) ((BoxCluster)box).rejectNotify();
    }
}
