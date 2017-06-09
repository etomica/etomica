/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator.mcmove;

import etomica.action.AtomActionTranslateBy;
import etomica.action.MoleculeChildAtomAction;
import etomica.box.Box;
import etomica.atom.IMolecule;
import etomica.potential.PotentialMaster;
import etomica.util.random.IRandom;
import etomica.simulation.Simulation;
import etomica.atom.MoleculeSource;
import etomica.atom.MoleculeSourceRandomMolecule;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.atom.iterator.MoleculeIterator;
import etomica.atom.iterator.MoleculeIteratorSinglet;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Standard Monte Carlo molecule-displacement trial move.
 *
 * @author David Kofke
 */
public class MCMoveMolecule extends MCMoveBoxStep implements MCMoveMolecular {
    
    private static final long serialVersionUID = 1L;
    protected final AtomIteratorArrayListSimple affectedAtomIterator = new AtomIteratorArrayListSimple();
    protected final MoleculeIteratorSinglet affectedMoleculeIterator = new MoleculeIteratorSinglet();
    protected final MeterPotentialEnergy energyMeter;
    protected double uOld;
    protected double uNew = Double.NaN;
    protected final IRandom random;
    protected Space space;

    protected final MoleculeChildAtomAction moveMoleculeAction;
    protected final Vector groupTranslationVector;
    protected MoleculeSource moleculeSource;
    protected IMolecule molecule;

    public MCMoveMolecule(Simulation sim, PotentialMaster potentialMaster,
                          Space _space) {
        this(potentialMaster, sim.getRandom(), _space, 1.0, 15.0);
    }
    
    public MCMoveMolecule(PotentialMaster potentialMaster, IRandom random,
                          Space _space, double stepSize,
                          double stepSizeMax) {
        super(potentialMaster);
        this.random = random;
        this.space = _space;
        moleculeSource = new MoleculeSourceRandomMolecule();
        ((MoleculeSourceRandomMolecule)moleculeSource).setRandomNumberGenerator(random);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        setStepSizeMax(stepSizeMax);
        setStepSizeMin(0.0);
        setStepSize(stepSize);
        perParticleFrequency = true;
        energyMeter.setIncludeLrc(false);

        AtomActionTranslateBy translator = new AtomActionTranslateBy(_space);
        groupTranslationVector = translator.getTranslationVector();
        moveMoleculeAction = new MoleculeChildAtomAction(translator);
        
        //set directive to exclude intramolecular contributions to the energy

        MoleculeSourceRandomMolecule randomMoleculeSource = new MoleculeSourceRandomMolecule();
        randomMoleculeSource.setRandomNumberGenerator(random);
        setMoleculeSource(randomMoleculeSource);
    }
    

    public boolean doTrial() {
        if(box.getMoleculeList().getMoleculeCount()==0) return false;
        
        
        molecule = moleculeSource.getMolecule();
        if (molecule == null) return false;
        energyMeter.setTarget(molecule);
        uOld = energyMeter.getDataAsScalar();
        if(Double.isInfinite(uOld)) {
            throw new RuntimeException("Started with overlap");
        }
        groupTranslationVector.setRandomCube(random);
        groupTranslationVector.TE(stepSize);
        moveMoleculeAction.actionPerformed(molecule);
        return true;
    }
    
    /**
     * Returns log of the ratio of the trial probabilities, ln(Tij/Tji) for the
     * states encountered before (i) and after (j) the most recent call to doTrial(). 
     * Tij is the probability that this move would generate state j from state i, and
     * Tji is the probability that a subsequent call to doTrial would return to state i
     * from state j.
     */
    public double getA() {return 1.0;}
    
    /**
     * Returns the log of the limiting-distribution probabilities of states, ln(Pj/Pi), 
     * for the states encountered before (i) and after (j) the most recent call to 
     * doTrial.
     */
    public double getB() {
        uNew = energyMeter.getDataAsScalar();
        //System.out.println("Translation uNew: "+uNew+ " uOld: "+uOld+" re "+energyMeter.getDataAsScalar());
        return -(uNew - uOld);
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
