/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.molecule;

import etomica.species.SpeciesAgentManager;
import etomica.box.BoxAgentManager;
import etomica.box.BoxAgentManager.BoxAgentSource;
import etomica.box.BoxCellManager;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.species.ISpecies;
import etomica.util.Arrays;

/**
 * @author taitan
 */
public abstract class PotentialMasterNbrMolecular extends PotentialMaster {
    private static final SpeciesAgentManager.AgentSource<PotentialArrayMolecular> speciesAgentSource = new SpeciesAgentManager.AgentSource<PotentialArrayMolecular>() {
        @Override
        public PotentialArrayMolecular makeAgent(ISpecies type) {
            return new PotentialArrayMolecular();
        }

        @Override
        public void releaseAgent(PotentialArrayMolecular agent, ISpecies type) {}
    };

    protected final SpeciesAgentManager<PotentialArrayMolecular> rangedAgentManager;
    protected final SpeciesAgentManager<PotentialArrayMolecular> speciesAgentManager;
    protected final Simulation simulation;
    protected IPotential[] allPotentials = new IPotential[0];
    protected BoxAgentSource<? extends BoxCellManager> boxAgentSource;
    protected BoxAgentManager<? extends BoxCellManager> boxAgentManager;

    protected PotentialMasterNbrMolecular(Simulation sim, BoxAgentSource<? extends BoxCellManager> boxAgentSource,
                                          BoxAgentManager<? extends BoxCellManager> boxAgentManager) {
        super();
        simulation = sim;
        this.boxAgentSource = boxAgentSource;
        this.boxAgentManager = boxAgentManager;
        rangedAgentManager = new SpeciesAgentManager<>(speciesAgentSource, sim);
        speciesAgentManager = new SpeciesAgentManager<>(speciesAgentSource, sim);
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
        for (int i = 0; i < species.length; i++) {
            addRangedPotential(potential, species[i]);
        }
        addRangedPotentialForSpecies(potential, species);
    }

    protected abstract void addRangedPotentialForSpecies(IPotentialMolecular subPotential, ISpecies[] species);

    protected void addRangedPotential(IPotentialMolecular potential, ISpecies species) {

        PotentialArrayMolecular potentialMoleculeSpecies = rangedAgentManager.getAgent(species);
        potentialMoleculeSpecies.addPotential(potential);
        boolean found = false;
        for (int i = 0; i < allPotentials.length; i++) {
            if (allPotentials[i] == potential) {
                found = true;
            }
        }
        if (!found) {
            allPotentials = Arrays.addObject(allPotentials, potential);
        }
    }

    public void removePotential(IPotentialMolecular potential) {
        super.removePotential(potential);
        if (potential.getRange() < Double.POSITIVE_INFINITY) {
            this.rangedAgentManager.getAgents().values().forEach(potentialArrayMolecular -> {
                potentialArrayMolecular.removePotential(potential);
            });
        } else if (potential instanceof PotentialGroup) {
            this.speciesAgentManager.getAgents().values().forEach(potentialArrayMolecular -> {
                potentialArrayMolecular.removePotential(potential);
            });
        }
        allPotentials = Arrays.removeObject(allPotentials, potential);
    }

    public PotentialArrayMolecular getRangedPotentials(ISpecies species) {
        return rangedAgentManager.getAgent(species);
    }

    public PotentialArrayMolecular getIntraPotentials(ISpecies atomType) {
        return speciesAgentManager.getAgent(atomType);
    }

    public final BoxAgentManager<? extends BoxCellManager> getCellAgentManager() {
        return boxAgentManager;
    }

    /**
     * Returns the simulation associated with this PotentialMaster
     */
    public Simulation getSimulation() {
        return simulation;
    }
}
