/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.paracetamol;

import etomica.action.AtomActionTranslateBy;
import etomica.action.MoleculeChildAtomAction;
import etomica.atom.AtomArrayList;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.molecule.IMolecule;
import etomica.molecule.MoleculeSource;
import etomica.molecule.MoleculeSourceRandomMolecule;
import etomica.potential.PotentialGroup;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;


/**
 * Standard Monte Carlo molecule-displacement trial move.  Two molecules are moved at a
 * time in such a way that the geometric center of the system is not changed.
 *
 * @author Nancy Cribbin
 */
public class MCMoveMoleculeCoupledDLPOLY extends MCMoveBoxStep {

    private static final long serialVersionUID = 1L;
    protected final MoleculeChildAtomAction moveMoleculeAction;
    protected final Vector groupTransVect;
    protected IMolecule molecule0, molecule1;
    protected final MeterPotentialEnergy energyMeter;
    protected MoleculeSource moleculeSource;
    protected double uOld, uNew;
    protected final IRandom random;
    protected final AtomIteratorArrayListSimple affectedMoleculeIterator;
    protected final AtomArrayList affectedMoleculeList;
    protected final AtomActionTranslateBy singleAction;
    protected PotentialGroup potential;
    
    public MCMoveMoleculeCoupledDLPOLY(PotentialMaster potentialMaster, IRandom nRandom,
                                       Space _space){
        super(potentialMaster);
        this.random = nRandom;
        moleculeSource = new MoleculeSourceRandomMolecule();
        ((MoleculeSourceRandomMolecule)moleculeSource).setRandomNumberGenerator(random);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        
        affectedMoleculeList = new AtomArrayList(2);
        affectedMoleculeIterator = new AtomIteratorArrayListSimple(affectedMoleculeList);
        
        singleAction = new AtomActionTranslateBy(_space);
        groupTransVect = singleAction.getTranslationVector();
        
        moveMoleculeAction = new MoleculeChildAtomAction(singleAction);
        
        perParticleFrequency = true;
        energyMeter.setIncludeLrc(false);
    }

    public void setBox(Box newBox) {
        super.setBox(newBox);
        moleculeSource.setBox(newBox);
        energyMeter.setBox(newBox);
    }
    
    public void setPotential(PotentialGroup newPotential){
        potential = newPotential;
    }
    
    public AtomIterator affectedAtoms() {
        affectedMoleculeList.clear();
        affectedMoleculeList.addAll(molecule0.getChildList());
        affectedMoleculeList.addAll(molecule1.getChildList());
        return affectedMoleculeIterator;
    }

    public double energyChange() {return uNew - uOld;}

    public void acceptNotify() {
        // I do believe nothing needs to happen here.
    }

    public boolean doTrial() {
//        System.out.println("doTrial MCMoveMoleculeCoupled called");
        
        molecule0 = moleculeSource.getMolecule();
        molecule1 = moleculeSource.getMolecule();
        if(molecule0==null || molecule1==null || molecule0==molecule1) return false;
        
        //make sure we don't double count the molecule0-molecule1 interaction
        uOld = energyMeter.getDataAsScalar();
        
        
        if(uOld > 1e10){
            throw new ConfigurationOverlapException(box);
        }
        
        groupTransVect.setRandomCube(random);
        groupTransVect.TE(stepSize);
        moveMoleculeAction.actionPerformed(molecule0);
        groupTransVect.TE(-1.0);
        moveMoleculeAction.actionPerformed(molecule1);
        
        uNew = energyMeter.getDataAsScalar();
        /*
         * Because we have uNew is infinity, and we don't want to have to 
         * worry about the system subtracting infinity from infinity, and
         * setting uNew equal to zero, and accepting the move.
         */
        return true;
    }

    public double getChi(double temperature) {
        return Math.exp(-(uNew - uOld) / temperature);
    }

    public void rejectNotify() {
        moveMoleculeAction.actionPerformed(molecule0);
        groupTransVect.TE(-1.0);
        moveMoleculeAction.actionPerformed(molecule1);
    }
    
}
