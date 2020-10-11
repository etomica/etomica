/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator.mcmove;

import etomica.action.AtomActionTranslateBy;
import etomica.action.MoleculeChildAtomAction;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.molecule.IMolecule;
import etomica.molecule.MoleculeSource;
import etomica.molecule.MoleculeSourceRandomMolecule;
import etomica.molecule.iterator.MoleculeIterator;
import etomica.molecule.iterator.MoleculeIteratorSinglet;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

/**
 * Standard Monte Carlo molecule-displacement trial move.
 *
 * @author David Kofke
 */
public class MCMoveMolecule extends MCMoveBoxStep implements MCMoveMolecular {
    
    protected final AtomIteratorArrayListSimple affectedAtomIterator = new AtomIteratorArrayListSimple();
    protected final MoleculeIteratorSinglet affectedMoleculeIterator = new MoleculeIteratorSinglet();
    protected final MeterPotentialEnergy energyMeter;
    protected double uOld;
    protected double uNew = Double.NaN;
    protected final IRandom random;
    protected final Space space;

    protected final MoleculeChildAtomAction moveMoleculeAction;
    protected final Vector groupTranslationVector;
    protected MoleculeSource moleculeSource;
    protected IMolecule molecule;

    public MCMoveMolecule(Simulation sim, PotentialMaster potentialMaster,
                          Space space) {
        this(potentialMaster, sim.getRandom(), space, 1.0, 15.0);
    }
    
    public MCMoveMolecule(PotentialMaster potentialMaster, IRandom random,
                          Space space, double stepSize,
                          double stepSizeMax) {
        super(potentialMaster);
        this.random = random;
        this.space = space;
        moleculeSource = new MoleculeSourceRandomMolecule();
        ((MoleculeSourceRandomMolecule)moleculeSource).setRandomNumberGenerator(random);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        setStepSizeMax(stepSizeMax);
        setStepSizeMin(0.0);
        setStepSize(stepSize);
        perParticleFrequency = true;
        energyMeter.setIncludeLrc(false);

        AtomActionTranslateBy translator = new AtomActionTranslateBy(space);
        groupTranslationVector = translator.getTranslationVector();
        moveMoleculeAction = new MoleculeChildAtomAction(translator);
        
        //set directive to exclude intramolecular contributions to the energy

        MoleculeSourceRandomMolecule randomMoleculeSource = new MoleculeSourceRandomMolecule();
        randomMoleculeSource.setRandomNumberGenerator(random);
        setMoleculeSource(randomMoleculeSource);
    }
    

    public boolean doTrial() {
        if(box.getMoleculeList().size()==0) return false;
        
        
        molecule = moleculeSource.getMolecule();
        if (molecule == null) return false;
        if(potential.isPotentialHard()) {
            uOld = 0.0;
        } else {
            energyMeter.setTarget(molecule);
            uOld = energyMeter.getDataAsScalar();
            if(Double.isInfinite(uOld)) {
                throw new RuntimeException("Started with overlap");
            }
        }
        groupTranslationVector.setRandomCube(random);
        groupTranslationVector.TE(stepSize);
        moveMoleculeAction.actionPerformed(molecule);
        return true;
    }

    public double getChi(double temperature) {
        energyMeter.setTarget(molecule);
        uNew = energyMeter.getDataAsScalar();
        return Math.exp(-(uNew - uOld) / temperature);
    }
    
    public double energyChange() {return uNew - uOld;}
    
    /**
     * Method called by IntegratorMC in the event that the most recent trial is accepted.
     */
    public void acceptNotify() {  /* do nothing */
    }
    
    public void rejectNotify() {
        groupTranslationVector.TE(-1);
        moveMoleculeAction.actionPerformed(molecule);
    }

    public void setBox(Box p) {
        super.setBox(p);
        energyMeter.setBox(p);
        moleculeSource.setBox(p);
    }

    public AtomIterator affectedAtoms() {
        affectedAtomIterator.setList(molecule.getChildList());
        return affectedAtomIterator;
    }
    
    public MoleculeIterator affectedMolecules(Box box) {
        affectedMoleculeIterator.setMolecule(molecule);
        return affectedMoleculeIterator;
    }
    
    /**
     * @return Returns the atomSource.
     */
    public MoleculeSource getAtomSource() {
        return moleculeSource;
    }
    /**
     * @param source The atomSource to set.
     */
    public void setMoleculeSource(MoleculeSource source) {
        moleculeSource = source;
    }
}
