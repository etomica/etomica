/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.mcmove;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.molecule.CenterOfMass;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.nbr.cell.NeighborCellManager;
import etomica.potential.compute.NeighborManager;
import etomica.potential.compute.NeighborManagerCell;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.RotationTensor3D;
import etomica.species.ISpecies;
import etomica.util.collections.IntArrayList;
import etomica.util.random.IRandom;
import etomica.virial.BoxCluster;
import etomica.virial.PotentialMoleculePairCell;

/**
 * Monte Carlo move for Mayer sampling that explore bond angle degrees of
 * freedom.  One bond angle is varied in each molecule.  Atoms on each side of
 * the bond are rotated around the middle atom of the bond angle with the
 * overall molecule COM fixed.
 *
 * Originally developed by Arpit for alkanes.
 */
public class MCMoveClusterAngleMulti extends MCMoveBoxStep {
    private final PotentialCompute potential;
    protected final IRandom random;
    protected final Space space;
    protected ISpecies species;

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
    protected boolean fixedCOM = true;
    //protected int start, stop;#delete
    protected boolean doLattice;
    protected int[] constraintMap;
    protected final MCMoveClusterAngle.AtomChooser atomChooser;
    protected NeighborManagerCell cellManager;
    public boolean skipW=false;
    protected final int numAtomsMoved;

    public MCMoveClusterAngleMulti(PotentialCompute potentialCompute, Space space, IntArrayList[] bonding, IRandom random, double stepSize) {
        this (potentialCompute,space,bonding,random,stepSize,new MCMoveClusterAngle.AtomChooserSimple(random),5);

    }
    public MCMoveClusterAngleMulti(PotentialCompute potentialCompute, Space space, IntArrayList[] bonding, IRandom random, double stepSize, MCMoveClusterAngle.AtomChooser atomChooser, int numAtomsMoved){
        super();
        this.potential = potentialCompute;
        this.space = space;
        this.random = random;
        this.bonding = bonding;
        this.stepSize = stepSize;
        this.atomChooser = atomChooser;
        modified = new int[bonding.length];
        setStepSizeMax(Math.PI);
        this.numAtomsMoved=numAtomsMoved;
    }
    public void setFixedCOM(boolean fixedCOM){
        this.fixedCOM=fixedCOM;

    }

    public void setDoLattice(boolean doLattice) {
        this.doLattice = doLattice;
    }

    public void setBox(Box p) {
        super.setBox(p);
        position = space.makeVectorArray(p.getMoleculeList().get(0).getChildList().size());
    }
    public void setCellManager(NeighborManagerCell cellManager) {
        this.cellManager = cellManager;
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
        if (!skipW && box instanceof BoxCluster) wOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        IMoleculeList moleculeList = box.getMoleculeList();
        iMolecule = random.nextInt(moleculeList.size());
        while (species != null && moleculeList.get(iMolecule).getType() != species) {
            iMolecule = random.nextInt(moleculeList.size());
        }

        IMolecule molecule = moleculeList.get(iMolecule);
        IAtomList atoms = molecule.getChildList();
        for(int i = 0 ;i<2 && false;i++){
            for (int j = 0;j<8; j++){
                int aa = 1+8*i+j ;
                int bb = aa - 1;
                if (j==0) bb=0;
               double r2=atoms.get(aa).getPosition().Mv1Squared(atoms.get(bb).getPosition());
                if(Math.abs(r2-1)>1e-6) throw new RuntimeException(aa+" "+bb+" "+r2+" ");
            }
        }
        for(int j = 0; j < molecule.getChildList().size(); j++) {
            position[j].E(atoms.get(j).getPosition());
        }
        if (molecule.getChildList().size() < 3) return false;
        modifiedIndex = 0;
        int[]ba=atomChooser.chooseAtoms(bonding);
        int a=ba[0];
        b=ba[1];


        modified[modifiedIndex]=b;
        ++modifiedIndex;

        modified[modifiedIndex] = a;
        ++modifiedIndex;
        Vector r = space.makeVector();
        r.Ev1Mv2(atoms.get(b).getPosition(), atoms.get(a).getPosition());
        Vector axis = space.makeVector();
        double dt=0;
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
        //System.out.println("axis " +axis+" dt "+dt);
        RotationTensor3D rotationTensor = new RotationTensor3D();
        rotationTensor.setRotationAxis(axis, dt);
        Vector shift = space.makeVector();
        //System.out.println("rotating "+a+" about "+b);
        transform(rotationTensor, a, atoms, shift);
        transformBondedAtoms(rotationTensor, a, atoms, shift);


        if(fixedCOM && (iMolecule==0 || (constraintMap != null && constraintMap[iMolecule] == 0))) {
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

        if (cellManager != null){
          cellManager.assignCellAll();
        }
        if (box instanceof BoxCluster) {
            ((BoxCluster)box).trialNotify();
            if (!skipW) wNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);

        }

        uNew = potential.computeAll(false);
for(int i = 0 ;i<2 && false;i++){
    for (int j = 0;j<8; j++){
        int aa = 1+8*i+j ;
        int bb = aa - 1;
        if (j==0) bb=0;
        double r2=atoms.get(aa).getPosition().Mv1Squared(atoms.get(bb).getPosition());
        if(Math.abs(r2-1)>1e-6) throw new RuntimeException(aa+" "+bb+" "+r2+" "+a+" "+b);
    }
}

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
                atoms.get(bonding[index].getInt(k)).getPosition().PE(shift);
                if(modifiedIndex-1<numAtomsMoved) {

                    double dt = 2 * stepSize * (random.nextDouble() - 0.5);
                    Vector r = space.makeVector();
                    r.Ev1Mv2(atoms.get(index).getPosition(), atoms.get(bonding[index].getInt(k)).getPosition());
                    Vector axis = space.makeVector();
                    axis.setRandomSphere(random);
                    Vector projection = space.makeVector();
                    projection.Ea1Tv1(axis.dot(r) / r.squared(), r);
                    axis.ME(projection);
                    axis.normalize();
                    rotationTensor3D.setRotationAxis(axis, dt);
                    b=index;
                    transform(rotationTensor3D, bonding[index].getInt(k), atoms, shift);
                }
                modified[modifiedIndex] = bonding[index].getInt(k);
                ++modifiedIndex;
                transformBondedAtoms(rotationTensor3D, bonding[index].getInt(k), atoms, shift);
            }
        }
    }

    protected void transform(RotationTensor3D rotationTensor3D, int index, IAtomList atoms, Vector shift) {
        //System.out.println(" rotating "+index);
        Vector r = space.makeVector();
        IAtom a = atoms.get(index);

        Vector p = a.getPosition();
        shift.ME( p);
        r.Ev1Mv2(p, atoms.get(b).getPosition());
        rotationTensor3D.transform(r);
        r.PE(atoms.get(b).getPosition());
        //System.out.println("r "+ r +" p "+p);
        p.E(r);
        shift.PE(p);
    }

    @Override
    public double getChi(double temperature) {
        if(IntegratorMC.dodebug && box.getIndex()==0){
            System.out.println("old "+wOld+" "+uOld+" new "+ wNew+" "+uNew);
        }
        return (wOld == 0 ? 1 : wNew / wOld) * Math.exp(-(uNew - uOld) / temperature);
    }

    @Override
    public void acceptNotify() {
        if (box instanceof BoxCluster) ((BoxCluster)box).acceptNotify();
    }

    @Override
    public void rejectNotify() {
        IAtomList atomList = box.getMoleculeList().get(iMolecule).getChildList();
        for(int j = 0; j < atomList.size(); j++) {
           atomList.get(j).getPosition().E(position[j]);

        }
        if(cellManager!= null){
            cellManager.assignCellAll();
        }
        if (box instanceof BoxCluster) ((BoxCluster)box).rejectNotify();
    }


}

