/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator.mcmove;

import etomica.action.AtomActionTranslateBy;
import etomica.action.MoleculeActionTranslateTo;
import etomica.action.MoleculeChildAtomAction;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.atom.iterator.AtomIteratorNull;
import etomica.box.Box;
import etomica.box.RandomPositionSource;
import etomica.box.RandomPositionSourceRectangular;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorMC;
import etomica.molecule.*;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

/**
 * Performs a trial that results in the exchange of a molecule from one box to another.
 * Primary use is as an elementary move in a Gibbs ensemble simulation
 *
 * @author David Kofke
 */
 
public class MCMoveMoleculeExchange extends MCMove {
    
    private static final long serialVersionUID = 2L;
    protected Box box1;
    protected Box box2;
    protected final IntegratorBox integrator1, integrator2;
    private final MeterPotentialEnergy energyMeter;
    private final AtomIteratorArrayListSimple affectedAtomIterator = new AtomIteratorArrayListSimple();
    private final MoleculeActionTranslateTo moleculeTranslator;
    private final MoleculeChildAtomAction moleculeReplacer;
    private final Vector translationVector;
    private final IRandom random;
    private MoleculeSource moleculeSource;
    protected RandomPositionSource positionSource;
    
    private transient IMolecule molecule;
    private transient Box iBox, dBox;
    private transient double uOld;
    private transient double uNew = Double.NaN;
    

    public MCMoveMoleculeExchange(PotentialMaster potentialMaster, IRandom random,
                                  Space space,
                                  IntegratorBox integrator1,
                                  IntegratorBox integrator2) {
        super(potentialMaster);
        this.random = random;
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        energyMeter.setIncludeLrc(true);
        moleculeReplacer = new MoleculeChildAtomAction(new AtomActionTranslateBy(space));
        moleculeTranslator = new MoleculeActionTranslateTo(space);
        translationVector = space.makeVector();
        setAtomPositionDefinition(new MoleculePositionCOM(space));
        this.integrator1 = integrator1;
        this.integrator2 = integrator2;
        box1 = integrator1.getBox();
        box2 = integrator2.getBox();
        moleculeSource = new MoleculeSourceRandomMolecule();
        ((MoleculeSourceRandomMolecule)moleculeSource).setRandomNumberGenerator(random);
        positionSource = new RandomPositionSourceRectangular(space, random);
    }
    
    public boolean doTrial() {
        if(random.nextInt(2) == 0) {
            iBox = box1;
            dBox = box2;
        }
        else {
            iBox = box2;
            dBox = box1;
        }
        if(dBox.getMoleculeList().getMoleculeCount() == 0) { //no molecules to delete; trial is over
            uNew = uOld = 0.0;
            return false;
        }

        moleculeSource.setBox(dBox);
        molecule = moleculeSource.getMolecule();  //select random molecule to delete

        if(potential.isPotentialHard()) {
            uOld = 0.0;
        } else {
            energyMeter.setBox(dBox);
            energyMeter.setTarget(molecule);
            uOld = energyMeter.getDataAsScalar();
        }
        dBox.removeMolecule(molecule);

        positionSource.setBox(iBox);
        moleculeTranslator.setDestination(positionSource.randomPosition());         //place at random in insertion box
        moleculeTranslator.actionPerformed(molecule);
        iBox.addMolecule(molecule);
        uNew = Double.NaN;
        return true;
    }//end of doTrial

    /**
     * Sets a new RandomPositionSource for this move to use.  By default, a
     * position source is used which assumes rectangular boundaries.
     */
    public void setPositionSource(RandomPositionSource newPositionSource) {
        positionSource = newPositionSource;
    }
    
    /**
     * Returns the RandomPositionSource used by this move.
     */
    public RandomPositionSource getPositionSource() {
        return positionSource;
    }

    /**
     * Sets the AtomSource this class uses to pick molecules to delete.
     */
    public void setMoleculeSource(MoleculeSource newMoleculeSource) {
        moleculeSource = newMoleculeSource;
    }
    
    /**
     * Returns the AtomSource this class uses to pick molecules to delete.
     */
    public MoleculeSource getMoleculeSource() {
        return moleculeSource;
    }

    public double getChi(double temperature) {
        energyMeter.setBox(iBox);
        energyMeter.setTarget(molecule);
        uNew = energyMeter.getDataAsScalar();
        double B = -(uNew - uOld);
        // assume both integrators have the same temperature
        //note that dSpecies.nMolecules has been decremented
        //and iSpecies.nMolecules has been incremented
        return Math.exp(B / temperature) * (dBox.getNMolecules(molecule.getType()) + 1) / dBox.getBoundary().volume()
               * iBox.getBoundary().volume()/iBox.getNMolecules(molecule.getType()); 
    }

    public void acceptNotify() {
        IntegratorBox iIntegrator = integrator1;
        IntegratorBox dIntegrator = integrator2;
        if (iIntegrator.getBox() == dBox) {
            iIntegrator = integrator2;
            dIntegrator = integrator1;
        }
        if (iIntegrator instanceof IntegratorMC) {
            ((IntegratorMC)iIntegrator).notifyEnergyChange(uNew);
        }
        else {
            //XXX grossly inefficient
            iIntegrator.reset();
        }
        if (dIntegrator instanceof IntegratorMC) {
            ((IntegratorMC)dIntegrator).notifyEnergyChange(-uOld);
        }
        else {
            //XXX grossly inefficient
            dIntegrator.reset();
        }
    }
    
    public void rejectNotify() {
        iBox.removeMolecule(molecule);
        translationVector.Ea1Tv1(-1,moleculeTranslator.getTranslationVector());
        ((AtomActionTranslateBy)moleculeReplacer.getAtomAction()).setTranslationVector(translationVector);
        moleculeReplacer.actionPerformed(molecule);
        dBox.addMolecule(molecule);
    }

    public final AtomIterator affectedAtoms(Box box) {
        if(this.box1 != box && this.box2 != box) return AtomIteratorNull.INSTANCE;
        affectedAtomIterator.setList(molecule.getChildList());
        return affectedAtomIterator;
    }
    
    public double energyChange(Box box) {
        if(box == iBox) return uNew;
        else if(box == dBox) return -uOld;
        else return 0.0;
    }

    
    /**
     * @return Returns the atomPositionDefinition.
     */
    public IMoleculePositionDefinition getAtomPositionDefinition() {
        return moleculeTranslator.getAtomPositionDefinition();
    }
    /**
     * @param atomPositionDefinition The atomPositionDefinition to set.
     */
    public void setAtomPositionDefinition(
            IMoleculePositionDefinition atomPositionDefinition) {
        moleculeTranslator.setAtomPositionDefinition(atomPositionDefinition);
    }
}//end of MCMoveMoleculeExchange
