/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator.mcmove;

import etomica.action.MoleculeActionTranslateTo;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.molecule.*;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.species.ISpecies;
import etomica.util.random.IRandom;

/**
 * Basic Monte Carlo move for semigrand-ensemble simulations.  Move consists
 * of selecting a molecule at random and changing its species identity.  More precisely,
 * the molecule is removed and another molecule of a different species replaces it.
 * An arbitrary number of species may be designated as subject to these exchange moves.
 * Acceptance is regulated by a set of fugacity fractions that are specified at design time.
 *
 * @author Jhumpa Adhikari
 * @author David Kofke
 */
public class MCMoveSemigrand extends MCMoveBox {
    
    private static final long serialVersionUID = 2L;
    private ISpecies[] speciesSet;
    private MoleculeArrayList[] reservoirs;
    private double[] fugacityFraction;
    private int nSpecies;
    private final AtomIteratorArrayListSimple affectedAtomIterator; 
    private final MeterPotentialEnergy energyMeter;
    private final MoleculeActionTranslateTo moleculeTranslator;
    private IMoleculePositionDefinition atomPositionDefinition;
    private final IRandom random;
    
    private transient IMolecule deleteMolecule, insertMolecule;
    private transient double uOld;
    private transient double uNew = Double.NaN;
    private transient int iInsert, iDelete;

    public MCMoveSemigrand(PotentialMaster potentialMaster, IRandom random,
                           Space _space) {
        super(potentialMaster);
        this.random = random;
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        affectedAtomIterator = new AtomIteratorArrayListSimple();
        perParticleFrequency = true;
        energyMeter.setIncludeLrc(true);
        moleculeTranslator = new MoleculeActionTranslateTo(_space);
        setAtomPositionDefinition(new MoleculePositionCOM(_space));
    }
    
    /**
     * Extends the superclass method to initialize the exchange-set species agents for the box.
     */
    public void setBox(Box p) {
        super.setBox(p);
        energyMeter.setBox(box);
    }//end setBox
    
    /**
     * Mutator method for the set of species that can participate in an exchange move.
     */
    public void setSpecies(ISpecies[] species) {
        nSpecies = species.length;
        if(nSpecies < 2) throw new IllegalArgumentException("Wrong size of species array in MCMoveSemigrand");
        speciesSet = new ISpecies[nSpecies];
        fugacityFraction = new double[nSpecies];
        reservoirs = new MoleculeArrayList[nSpecies];
        for(int i=0; i<nSpecies; i++) {
            speciesSet[i] = species[i];
            fugacityFraction[i] = 1.0/nSpecies;
            reservoirs[i] = new MoleculeArrayList();
        }
    }
    
    /**
     * Accessor method for the set of species that can participate in an exchange move.
     */
    public ISpecies[] getSpecies() {return speciesSet;}
    
    /**
     * Specifies the fugacity fractions for the set of species that can participate in
     * an exchange move.  The given array must have the same dimension as the array of
     * species that was previously set in a call to setSpecies.  If the given set of "fractions"
     * does not sum to unity, the values will be normalized (e.g., sending the set {1.0, 1.0} 
     * leads to fugacity fractions of {0.5, 0.5}).
     */
    public void setFugacityFraction(double[] f) {
        if(f.length != nSpecies || speciesSet == null) 
            throw new IllegalArgumentException("Wrong size of fugacity-fraction array in MCMoveSemigrand");
            
        double sum = 0.0;
        for(int i=0; i<nSpecies; i++) {
            fugacityFraction[i] = f[i]; 
            if(f[i] < 0.0) throw new IllegalArgumentException("Negative fugacity-fraction MCMoveSemigrand");
            sum += f[i];
        }
        for(int i=0; i<nSpecies; i++) {fugacityFraction[i] /= sum;}//normalize to unity
    }

    public double getFugacityFraction(int i) {
        if(i < 0 || i >= nSpecies) 
            throw new IllegalArgumentException("Illegal fugacity-fraction index in MCMoveSemigrand");
        return fugacityFraction[i];
    }

    /**
     * Accessor method for the set of fugacity fractions.
     */
    public double[] getFugacityFraction() {return fugacityFraction;}
    
    public boolean doTrial() {
        //select species for deletion
        iDelete = random.nextInt(nSpecies);//System.out.println("Random no. :"+randomNo);
        if(box.getNMolecules(speciesSet[iDelete]) == 0) {
            uNew = uOld = 0.0;
            return false;
        }

        //select species for insertion
        iInsert = iDelete;
        if(nSpecies == 2) iInsert = 1 - iDelete;
        else while(iInsert == iDelete) {iInsert = random.nextInt(nSpecies);}
  
        IMoleculeList moleculeList = box.getMoleculeList(speciesSet[iDelete]);
        deleteMolecule = moleculeList.getMolecule(random.nextInt(moleculeList.getMoleculeCount()));
        energyMeter.setTarget(deleteMolecule);
        uOld = energyMeter.getDataAsScalar();
        box.removeMolecule(deleteMolecule);
        
        int size = reservoirs[iInsert].getMoleculeCount();
        if(size>0) {
            insertMolecule = reservoirs[iInsert].remove(size-1);
            box.addMolecule(insertMolecule);
        }
        else {
            insertMolecule = speciesSet[iInsert].makeMolecule();
            box.addMolecule(insertMolecule);
        }
        moleculeTranslator.setDestination(atomPositionDefinition.position(deleteMolecule));
        moleculeTranslator.actionPerformed(insertMolecule);
        //in general, should also randomize orintation and internal coordinates
        uNew = Double.NaN;
        return true;
    }//end of doTrial

    public double getChi(double temperature) {
        energyMeter.setTarget(insertMolecule);
        uNew = energyMeter.getDataAsScalar();
        double B = Math.exp(-(uNew - uOld));
        return (double) (box.getNMolecules(speciesSet[iDelete]) + 1) / (double) box.getNMolecules(speciesSet[iInsert])
                * (fugacityFraction[iInsert] / fugacityFraction[iDelete]) * B;
    }

    public void acceptNotify() {
        //put deleted molecule in reservoir
        reservoirs[iDelete].add(deleteMolecule);
    }

    public void rejectNotify() {
        //put deleted molecule back into box
        box.addMolecule(deleteMolecule);
        //remove inserted molecule and put in reservoir
        box.removeMolecule(insertMolecule);
        reservoirs[iInsert].add(insertMolecule);
    }
    
    

    public double energyChange() {return uNew - uOld;}
    
    public final AtomIterator affectedAtoms() {
        affectedAtomIterator.setList(insertMolecule.getChildList());
        return affectedAtomIterator;
    }

    /**
     * @return Returns the positionDefinition.
     */
    public IMoleculePositionDefinition geAtomPositionDefinition() {
        return atomPositionDefinition;
    }

    /**
     * @param positionDefinition The positionDefinition to set.
     */
    public void setAtomPositionDefinition(IMoleculePositionDefinition positionDefinition) {
        this.atomPositionDefinition = positionDefinition;
        moleculeTranslator.setAtomPositionDefinition(positionDefinition);
    }

}
