/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.IAtomList;
import etomica.atom.iterator.AtomIterator;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.space3d.Vector3D;
import etomica.util.random.IRandom;

/**
 * An MC move for cluster simulations which performs reptation moves on molecules.
 * One of the atoms (the first or last atom) is moved to a random position
 * 1 bondlength away from its current position and the other Atoms are moved
 * into the position of the adjacent Atom.
 *
 * @author Andrew Schultz
 */
public class MCMoveClusterReptateMulti extends MCMoveBox {

    private static final long serialVersionUID = 2L;
    private final MeterPotentialEnergy energyMeter;
    protected final IRandom random;

    public MCMoveClusterReptateMulti(Simulation sim, PotentialMaster potentialMaster, int nAtoms) {
    	this(potentialMaster, sim.getRandom(), nAtoms);
        setBondLength(1.0);
    }
    
    public MCMoveClusterReptateMulti(PotentialMaster potentialMaster, IRandom random, int nAtoms) {
        super(potentialMaster);
        this.random = random;
        this.nAtoms = nAtoms;
        selectedMolecules = new IMolecule[nAtoms];
        oldPositions = new Vector3D[nAtoms];
        for (int i=0; i<nAtoms; i++) {
            oldPositions[i] = new Vector3D();
        }
        forward = new boolean[nAtoms];
        energyMeter = new MeterPotentialEnergy(potential);
        work1 = new Vector3D();
    }

    public void setBox(Box p) {
        super.setBox(p);
        energyMeter.setBox(p);
    }
    
    //note that total energy is calculated
    public boolean doTrial() {
        if (selectedMolecules[0] == null) selectMolecules();
//        System.out.println("old energy");
//        Potential2HardSpherical.foo = true;
        uOld = energyMeter.getDataAsScalar();
        if (Double.isInfinite(uOld)) {
            uOld = energyMeter.getDataAsScalar();
            throw new IllegalStateException("U can't be infinite before the move");
        }
//        Potential2HardSpherical.foo = false;
//        System.out.println("old energy done");
        wOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        for(int i=0; i<selectedMolecules.length; i++) {
            forward[i] = random.nextInt(2) == 0;
            IAtomList childList = selectedMolecules[i].getChildList();
            int numChildren = childList.size();
            for (int k=0; k<numChildren; k++) {
//                System.out.println(i+" before "+k+" "+((AtomLeaf)childList.get(k)).coord.position());
                if (k > 0) {
                    work1.E(childList.get(k).getPosition());
                    work1.ME(childList.get(k-1).getPosition());
                    double d = Math.sqrt(work1.squared());
//                    System.out.println("distance "+d);
                    if (Math.abs(d - bondLength)/bondLength > 0.0000001) {
                        throw new IllegalStateException("wiggle "+i+" "+k+" bond length should be close to "+bondLength+" ("+d+")");
                    }
                }
            }
            if (forward[i]) {
                Vector position = childList.get(numChildren-1).getPosition();
                oldPositions[i].E(position);
                for (int j=numChildren-1; j>0; j--) {
                    Vector position2 = childList.get(j-1).getPosition();
                    position.E(position2);
                    position = position2;
                }
                work1.setRandomSphere(random);
                work1.TE(bondLength);
                childList.get(0).getPosition().PE(work1);
            }
            else {
                Vector position = childList.get(0).getPosition();
                oldPositions[i].E(position);
                for (int j=0; j<numChildren-1; j++) {
                    Vector position2 = childList.get(j+1).getPosition();
                    position.E(position2);
                    position = position2;
                }
                work1.setRandomSphere(random);
                work1.TE(bondLength);
                childList.get(numChildren-1).getPosition().PE(work1);
            }

            for (int k=0; k<numChildren; k++) {
//                System.out.println(i+" after "+k+" "+((AtomLeaf)childList.get(k)).coord.position());
                if (k > 0) {
                    work1.E(childList.get(k).getPosition());
                    work1.ME(childList.get(k-1).getPosition());
                    double d = Math.sqrt(work1.squared());
//                    System.out.println("distance "+d);
                    if (Math.abs(d - bondLength)/bondLength > 0.0000001) {
                        throw new IllegalStateException("wiggle "+i+" "+k+" bond length should be close to "+bondLength+" ("+d+")");
                    }
                }
            }
        }
        ((BoxCluster)box).trialNotify();
        wNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
//        System.out.println("now energy");
//        Potential2HardSpherical.foo = true;
        uNew = energyMeter.getDataAsScalar();
//        Potential2HardSpherical.foo = false;
//        System.out.println(uOld+" => "+uNew+"   "+wOld+" => "+wNew);
        return true;
    }
    
    public void setBondLength(double b) {
        bondLength = b;
    }
	
    protected void selectMolecules() {
        IMoleculeList moleculeList = box.getMoleculeList();
        if (moleculeList.getMoleculeCount() != nAtoms+1) throw new IllegalStateException("move should work on number of molecules in box - 1");
        //skip the first one
        for (int i=1; i<moleculeList.getMoleculeCount(); i++) {
            selectedMolecules[i++] = moleculeList.getMolecule(i);
        }
    }
	
    public void rejectNotify() {
        for(int i=0; i<selectedMolecules.length; i++) {
            IAtomList childList = selectedMolecules[i].getChildList();
            int numChildren = childList.size();
            if (!forward[i]) {
                Vector position = childList.get(numChildren-1).getPosition();
                for (int j=numChildren-1; j>0; j--) {
                    Vector position2 = childList.get(j-1).getPosition();
                    position.E(position2);
                    position = position2;
                }
                childList.get(0).getPosition().E(oldPositions[i]);
            }
            else {
                Vector position = childList.get(0).getPosition();
                for (int j=0; j<numChildren-1; j++) {
                    Vector position2 = childList.get(j+1).getPosition();
                    position.E(position2);
                    position = position2;
                }
                childList.get(numChildren-1).getPosition().E(oldPositions[i]);
            }
//            System.out.println("rejected");
        }
        ((BoxCluster)box).rejectNotify();
    }

    public void acceptNotify() {
        ((BoxCluster)box).acceptNotify();
    }
    
    public double getB() {
        return -(uNew - uOld);
    }

    public double getChi(double temperature) {
        return ((wOld == 0.0) ? 1 : wNew / wOld) * Math.exp(-(uNew - uOld) / temperature);
    }
    
    public double energyChange() {
        return uNew - uOld;
    }

    public AtomIterator affectedAtoms() {
        return null;
    }
	
    private final int nAtoms;
    private final IMolecule[] selectedMolecules;
    private double bondLength;
    private final Vector3D work1;
    private final Vector3D[] oldPositions;
    private final boolean[] forward;
    private double wOld, wNew;
    private double uOld, uNew;
}
