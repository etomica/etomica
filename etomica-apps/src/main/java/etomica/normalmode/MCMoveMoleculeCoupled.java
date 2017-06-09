/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.AtomActionTranslateBy;
import etomica.action.MoleculeChildAtomAction;
import etomica.atom.AtomArrayList;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.molecule.*;
import etomica.nbr.list.molecule.PotentialMasterListMolecular;
import etomica.potential.IPotentialMolecular;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

/**
 * Standard Monte Carlo molecule-displacement trial move.  Two molecules are moved at a
 * time in such a way that the geometric center of the system is not changed.
 *
 * @author Nancy Cribbin
 * 	       Tai Boon Tan
 */
public class MCMoveMoleculeCoupled extends MCMoveBoxStep {

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
    protected final MoleculePair pair;
    protected IPotentialMolecular pairPotential;
    protected boolean doExcludeNonNeighbors, doIncludePair;
    
    public MCMoveMoleculeCoupled(PotentialMaster potentialMaster, IRandom nRandom,
                                 Space _space){
        super(potentialMaster);
        this.random = nRandom;
        moleculeSource = new MoleculeSourceRandomMolecule();
        ((MoleculeSourceRandomMolecule)moleculeSource).setRandomNumberGenerator(random);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        
        affectedMoleculeList = new AtomArrayList();
        affectedMoleculeIterator = new AtomIteratorArrayListSimple(affectedMoleculeList);
        
        singleAction = new AtomActionTranslateBy(_space);
        groupTransVect = singleAction.getTranslationVector();
        
        moveMoleculeAction = new MoleculeChildAtomAction(singleAction);
        
        pair = new MoleculePair();
        
        perParticleFrequency = true;
        //energyMeter.setIncludeLrc(false);
    }

    public void setBox(Box newBox) {
        super.setBox(newBox);
        moleculeSource.setBox(newBox);
        energyMeter.setBox(newBox);
    }
    
    public void setPotential(IPotentialMolecular newPotential){
        pairPotential = newPotential;
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
        pair.atom0 = molecule0;
        pair.atom1 = molecule1;
        energyMeter.setTarget(molecule0);
        uOld = energyMeter.getDataAsScalar();
        energyMeter.setTarget(molecule1);
        uOld += energyMeter.getDataAsScalar();
        
        doIncludePair = pairPotential != null;
        
        if (doIncludePair && doExcludeNonNeighbors && potential instanceof PotentialMasterListMolecular) {
            doIncludePair = false;
            IMoleculeList[] list0 = ((PotentialMasterListMolecular)potential).getNeighborManager(box).getDownList(molecule0);
            for (int i=0; i<list0.length; i++) {
                if (((MoleculeArrayList)list0[i]).indexOf(molecule1)>-1) {
                    doIncludePair = true;
                    break;
                }
            }
            if (!doIncludePair) {
                list0 = ((PotentialMasterListMolecular)potential).getNeighborManager(box).getUpList(molecule0);
                for (int i=0; i<list0.length; i++) {
                    if (((MoleculeArrayList)list0[i]).indexOf(molecule1)>-1) {
                        doIncludePair = true;
                        break;
                    }
                }
            }
        }
        
        if(doIncludePair){
        	uOld -= pairPotential.energy(pair);
        }
             
        if(uOld > 1e10){
            throw new ConfigurationOverlapException(box);
        }
        
        groupTransVect.setRandomCube(random);
        groupTransVect.TE(stepSize);
        moveMoleculeAction.actionPerformed(molecule0);
        groupTransVect.TE(-1.0);
        moveMoleculeAction.actionPerformed(molecule1);
    
        return true;
    }

    public double getA() {
        return 1.0;
    }

    public double getB() {
        uNew = energyMeter.getDataAsScalar();
        energyMeter.setTarget(molecule0);
        uNew += energyMeter.getDataAsScalar();
        if(!Double.isInfinite(uNew) && doIncludePair){
            uNew -= pairPotential.energy(pair);
        }
        
        return -(uNew - uOld);
    }

    public void rejectNotify() {
        moveMoleculeAction.actionPerformed(molecule0);
        groupTransVect.TE(-1.0);
        moveMoleculeAction.actionPerformed(molecule1);
    }
    
    /**
     * Configures the move to not explicitly calculate the potential for atom
     * pairs that are not neighbors (as determined by the potentialMaster).
     * 
     * Setting this has an effect only if the potentialMaster is an instance of
     * PotentialMasterList
     */
    public void setDoExcludeNonNeighbors(boolean newDoExcludeNonNeighbors) {
        doExcludeNonNeighbors = newDoExcludeNonNeighbors;
    }

    /**
     * Returns true if the move does not explicitly calculate the potential for
     * atom pairs that are not neighbors (as determined by the potentialMaster).
     */
    public boolean getDoExcludeNonNeighbors() {
        return doExcludeNonNeighbors;
    }
    
}
