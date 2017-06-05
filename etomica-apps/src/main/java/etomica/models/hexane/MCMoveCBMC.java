/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.hexane;

import etomica.api.*;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.atom.AtomArrayList;
import etomica.atom.MoleculeSource;
import etomica.atom.MoleculeSourceRandomMolecule;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.potential.PotentialMaster;
import etomica.space.Vector;
import etomica.space.Space;
import etomica.util.Constants;
import etomica.util.random.IRandom;

public abstract class MCMoveCBMC extends MCMoveBox {

    public MCMoveCBMC(PotentialMaster potentialMaster, IRandom random,
                      Space _space, IntegratorMC integrator, Box p, int maxAtomsPerMolecule,
                      int NTrial) {
        super(potentialMaster);
        this.random = random;

        setNumberOfTrials(NTrial);

        beta = 1.0 / integrator.getTemperature() / Constants.BOLTZMANN_K;
        atomList = new AtomArrayList(maxAtomsPerMolecule);
        affectedAtomIterator = new AtomIteratorArrayListSimple();

        externalMeter = new MeterPotentialEnergy(potentialMaster);

        box = p;

        moleculeSource = new MoleculeSourceRandomMolecule();
        moleculeSource.setBox(box);
        ((MoleculeSourceRandomMolecule) moleculeSource).setRandomNumberGenerator(random);
        setMoleculeSource(moleculeSource);

        positionOld = new Vector[maxAtomsPerMolecule];
        for (int i = 0; i < maxAtomsPerMolecule; i++) {
            positionOld[i] = _space.makeVector();
        }
    }

    /**
     * Sets the AtomSource used to select molecules acted on by MC trials.
     */
    public void setMoleculeSource(MoleculeSource newMoleculeSource) {
        moleculeSource = newMoleculeSource;
    }

    /**
     * Returns the AtomSource used to select Atoms acted on by MC trials.
     */
    public MoleculeSource getMoleculeSource() {
        return moleculeSource;
    }

    public abstract double energyChange();

    public void setBox(Box p) {
        super.setBox(p);
        externalMeter.setBox(p);
    }

    public void acceptNotify() {
//         System.out.println("ACCEPTED A WHOLE MOVE!!!!!!!!!!!!!!!!!!!!!!");
    }

    public boolean doTrial() {
//        System.out.println("doTrial() CBMC called"); 

        // pick a molecule & get its childlist
        atom = moleculeSource.getMolecule();
        if (atom == null){
            return false;
        }
        affectedAtomIterator.setList(atom.getChildList());

        // we assume that that atoms that make the molecule are children of the
        // molecule.
        atomList = atom.getChildList();
        chainlength = atomList.getAtomCount();

        // store the old locations of every atom in the molecule in positionOld.
        for (int i = 0; i < chainlength; i++) {
            positionOld[i].E(atomList.getAtom(i).getPosition());
        }

        return calcRosenbluthFactors(); // this means we were able to propose a move.
    }

    public boolean doTrial(IMolecule atom){
        // who calls me?
        if (atom == null)
            return false;
        affectedAtomIterator.setList(atom.getChildList());

        // we assume that that atoms that make the molecule are children of the
        // molecule.
        atomList = atom.getChildList();
        chainlength = atomList.getAtomCount();

        // store the old locations of every atom in the molecule in positionOld.
        for (int i = 0; i < chainlength; i++) {
            positionOld[i].E(atomList.getAtom(i).getPosition());
        }

        calcRosenbluthFactors();
        return true; // this means we were able to propose a move.
    }
    
    public double getA() {
        return wNew / wOld;
    }

    public double getB() {
        return 0.0;
    }

    protected abstract boolean calcRosenbluthFactors();

    public void setNumberOfTrials(int n) {
        numTrial = n;
    }

    protected void setChainlength(int n) {
        chainlength = n;
    }

    public int getNumberOfTrials() {
        return numTrial;
    }

    public int getChainlength() {
        return chainlength;
    }

    public void rejectNotify() {
        for (int i = 0; i < chainlength; i++) {
             atomList.getAtom(i).getPosition().E(positionOld[i]);
        }
//        System.out.println("MCMoveCBMC rejects another!!");
    }

    public AtomIterator affectedAtoms() {
        return affectedAtomIterator;
    }

    
    private static final long serialVersionUID = 1L;

    protected final MeterPotentialEnergy externalMeter;

    protected double wNew; // Rosenbluth factor of new configuration

    protected double wOld; // Rosenbluth factor of old configuration

    protected int chainlength; // the number of atoms in a molecule; some
                                // juggling may be necessary to make this work
                                // if we want to make the chains longer....

    protected double beta;

    protected IMolecule atom;

    protected double uOld;

    protected double uNew = Double.NaN;

    protected Vector[] positionOld; // Used to store the position of the
                                        // molecule before mofing it.

    protected IAtomList atomList;

    protected int numTrial;

    protected MoleculeSource moleculeSource;

    private AtomIteratorArrayListSimple affectedAtomIterator;

    protected final IRandom random;

}
