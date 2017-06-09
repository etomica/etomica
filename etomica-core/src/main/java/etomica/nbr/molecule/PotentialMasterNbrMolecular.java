/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.molecule;

import etomica.atom.SpeciesAgentManager;
import etomica.box.BoxAgentManager;
import etomica.box.BoxAgentManager.BoxAgentSource;
import etomica.box.BoxCellManager;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.species.ISpecies;
import etomica.util.Arrays;

/**
 * @author taitan
 *
 */
public abstract class PotentialMasterNbrMolecular extends PotentialMaster implements SpeciesAgentManager.AgentSource {

    /**
	 * 
	 */

	protected PotentialMasterNbrMolecular(Simulation sim, BoxAgentSource<? extends BoxCellManager> boxAgentSource,
                                          BoxAgentManager<? extends BoxCellManager> boxAgentManager) {
        super();
        simulation = sim;
        this.boxAgentSource = boxAgentSource;
        this.boxAgentManager = boxAgentManager;
        rangedAgentManager = new SpeciesAgentManager(this);
        speciesAgentManager = new SpeciesAgentManager(this);

        rangedAgentManager.init(sim);
        speciesAgentManager.init(sim);
        
        rangedPotentialIterator = rangedAgentManager.makeIterator();
        speciesPotentialIterator = speciesAgentManager.makeIterator();
        boxAgentManager.setSimulation(sim);
    }
    
    public PotentialGroup makePotentialGroup(int nBody) {
    	throw new RuntimeException("<PotentialMasterNbrMolecular> The class does not recognize PotentialGroup!!!");
    }
    
    public void addPotential(IPotentialMolecular potential, ISpecies[] species) {
        if ((potential instanceof PotentialGroup)) {
            System.err.println("You gave me a concrete molecule potential and I'm very confused now.  I'll pretend like that's OK but don't hold your breath.");
            throw new RuntimeException("<PotentialMasterNbrMolecular> NO POTENTIAL GROUP in this class!!");
        }
        super.addPotential(potential, species);

        if (potential.getRange() == Double.POSITIVE_INFINITY && potential.nBody() > 1) {
            // -- should only happen for 0 or 1-body potentials, which should be fine
            throw new RuntimeException("<PotentialMasterNbrMolecular> infinite-ranged 2-body potential!!!");
        }
        for (int i=0; i<species.length; i++) {
            addRangedPotential(potential,species[i]);
        }
        addRangedPotentialForSpecies(potential, species);
    }

    

    protected abstract void addRangedPotentialForSpecies(IPotentialMolecular subPotential, ISpecies[] species);
    
    protected void addRangedPotential(IPotentialMolecular potential, ISpecies species) {
        
        PotentialArrayMolecular potentialMoleculeSpecies = (PotentialArrayMolecular)rangedAgentManager.getAgent(species);
        potentialMoleculeSpecies.addPotential(potential);
        boolean found = false;
        for (int i=0; i<allPotentials.length; i++) {
            if (allPotentials[i] == potential) {
                found = true;
            }
        }
        if (!found) {
            allPotentials = (IPotential[])etomica.util.Arrays.addObject(allPotentials, potential);
        }
    }
    
    public void removePotential(IPotentialMolecular potential) {
        super.removePotential(potential);
        if (potential.getRange() < Double.POSITIVE_INFINITY) {
            rangedPotentialIterator.reset();
            while (rangedPotentialIterator.hasNext()) {
                ((PotentialArrayMolecular)rangedPotentialIterator.next()).removePotential(potential);
            }
        }
        else if (potential instanceof PotentialGroup) {
            speciesPotentialIterator.reset();
            while (speciesPotentialIterator.hasNext()) {
                ((PotentialArrayMolecular)speciesPotentialIterator.next()).removePotential(potential);
            }
        }
        allPotentials = (IPotential[])Arrays.removeObject(allPotentials,potential);
    }
    
    public PotentialArrayMolecular getRangedPotentials(ISpecies species) {
        return (PotentialArrayMolecular)rangedAgentManager.getAgent(species);
    }

    public PotentialArrayMolecular getIntraPotentials(ISpecies atomType) {
        return (PotentialArrayMolecular)speciesAgentManager.getAgent(atomType);
    }
    
    public final BoxAgentManager<? extends BoxCellManager> getCellAgentManager() {
        return boxAgentManager;
    }
    
    public Class getSpeciesAgentClass() {
        return PotentialArrayMolecular.class;
    }
    
    public Object makeAgent(ISpecies species) {
        return new PotentialArrayMolecular();
    }
    
    public void releaseAgent(Object agent, ISpecies type) {
    }

    /**
     * Returns the simulation associated with this PotentialMaster
     */
    public Simulation getSimulation() {
        return simulation;
    }
    
	protected SpeciesAgentManager.AgentIterator rangedPotentialIterator;
    protected SpeciesAgentManager.AgentIterator speciesPotentialIterator;
    protected final SpeciesAgentManager rangedAgentManager;
    protected final SpeciesAgentManager speciesAgentManager;
    protected IPotential[] allPotentials = new IPotential[0];
    protected BoxAgentSource<? extends BoxCellManager> boxAgentSource;
    protected final Simulation simulation;
    protected BoxAgentManager<? extends BoxCellManager> boxAgentManager;
}
