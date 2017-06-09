/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator.mcmove;

import etomica.action.MoleculeActionTranslateTo;
import etomica.api.ISpecies;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.atom.iterator.AtomIteratorNull;
import etomica.box.Box;
import etomica.box.RandomPositionSource;
import etomica.box.RandomPositionSourceRectangular;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeArrayList;
import etomica.molecule.MoleculeSetSinglet;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.util.Debug;
import etomica.util.random.IRandom;

/**
 * Elementary Monte Carlo move in which a molecule of a specified species is
 * inserted into or removed from a box.
 *
 * @author David Kofke
 */
public class MCMoveInsertDelete extends MCMoveBox {

    //chemical potential
    protected double mu;
    
    //directive must specify "BOTH" to get energy with all atom pairs
    protected final MeterPotentialEnergy energyMeter;
	protected ISpecies species;
	protected final AtomIteratorArrayListSimple affectedAtomIterator = new AtomIteratorArrayListSimple();
	protected IMolecule testMolecule;
	protected double uOld;
	protected double uNew = Double.NaN;
	protected boolean insert;
	protected final MoleculeArrayList reservoir;
    protected final MoleculeActionTranslateTo atomTranslator;
    protected IMoleculeList moleculeList;
    protected IRandom random;
    protected RandomPositionSource positionSource;

    public MCMoveInsertDelete(PotentialMaster potentialMaster, IRandom random,
                              Space _space) {
        super(potentialMaster);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        setMu(0.0);
        energyMeter.setIncludeLrc(true);
        atomTranslator = new MoleculeActionTranslateTo(_space);
        reservoir = new MoleculeArrayList();
        this.random = random;
        perParticleFrequency = true;
        positionSource = new RandomPositionSourceRectangular(_space, random);
    }
    
//perhaps should have a way to ensure that two instances of this class aren't assigned the same species
    public void setSpecies(ISpecies s) {
        species = s;
        if(box != null) {
            moleculeList = box.getMoleculeList(species);
        }
    }
    public ISpecies getSpecies() {return species;}
    
    public void setBox(Box p) {
        super.setBox(p);
        energyMeter.setBox(box);
        if(species != null) {
            moleculeList = box.getMoleculeList(species);
        }
        positionSource.setBox(box);
    }
    
    /**
     * Sets a new RandomPositionSource for this move to use.  By default, a
     * position source is used which assumes rectangular boundaries.
     */
    public void setPositionSource(RandomPositionSource newPositionSource) {
        positionSource = newPositionSource;
        if (box != null) {
            positionSource.setBox(box);
        }
    }
    
    /**
     * Returns the RandomPositionSource used by this move.
     */
    public RandomPositionSource getPositionSource() {
        return positionSource;
    }

    /**
     * Chooses and performs with equal probability an elementary molecule insertion
     * or deletion.
     */
    public boolean doTrial() {
        insert = (random.nextInt(2) == 0);
        if(insert) {
            uOld = 0.0;
            
            if(!reservoir.isEmpty()) testMolecule = reservoir.remove(reservoir.getMoleculeCount()-1);
            else testMolecule = species.makeMolecule();

            atomTranslator.setDestination(positionSource.randomPosition());
            atomTranslator.actionPerformed(testMolecule);
            box.addMolecule(testMolecule);
        } else {//delete
            if(box.getNMolecules(species) == 0) {
                testMolecule = null;
                return false;
            }
            testMolecule = moleculeList.getMolecule(random.nextInt(moleculeList.getMoleculeCount()));
            //delete molecule only upon accepting trial
            energyMeter.setTarget(testMolecule);
            uOld = energyMeter.getDataAsScalar();
            uNew = 0.0;
        } 
        uNew = Double.NaN;
        return true;
    }//end of doTrial
    
    public double getA() {//note that moleculeCount() gives the number of molecules after the trial is attempted
        int numMolecules = box.getNMolecules(species);
        double a = box.getBoundary().volume()/numMolecules;
        return insert ? a : 1.0/a;
    }
    
    public double getB() {
        if(insert) {
            energyMeter.setTarget(testMolecule);
            uNew = energyMeter.getDataAsScalar();
        }
        else {
            uNew = 0;
        }
        double b = uOld - uNew;
        if (insert) b += mu;
        else b -= mu;
        return b;
    }
    
    public void acceptNotify() {
        if (insert && Debug.ON && Debug.DEBUG_NOW && Debug.allAtoms(new MoleculeSetSinglet(testMolecule))) System.out.println("accepted insertion of "+testMolecule);
        if(!insert) {
            // accepted deletion - remove from box and add to reservoir
            box.removeMolecule(testMolecule);
            reservoir.add(testMolecule);
        }
    }
    
    public void rejectNotify() {
        if(insert) {
            // rejected insertion - remove from box and return to reservoir
            box.removeMolecule(testMolecule);
            reservoir.add(testMolecule);
            // test molecule is no longer in the simulation and should not be 
            // returned by affectedAtoms
            testMolecule = null;
        }
    }
    
    public double energyChange() {return uNew - uOld;}

    /**
     * Returns an iterator giving molecule that is being added or deleted 
     * in the current or most recent trial.
     */
    public AtomIterator affectedAtoms() {
        if(testMolecule == null) return AtomIteratorNull.INSTANCE;
        affectedAtomIterator.setList(testMolecule.getChildList());
        return affectedAtomIterator;
    }

    public boolean lastMoveInsert() {
        return insert;
    }

    /**
     * Mutator method for the chemical potential of the insertion/deletion species.
     */
    public void setMu(double mu) {this.mu = mu;}
    
    /**
     * Accessor method for the chemical potential of th insertion/deletion species.
     */
    public final double getMu() {return mu;}
    
    /**
     * Indicates that chemical potential has dimensions of energy.
     */
    public final etomica.units.Dimension getMuDimension() {return etomica.units.Energy.DIMENSION;}
  
}//end of MCMoveInsertDelete
