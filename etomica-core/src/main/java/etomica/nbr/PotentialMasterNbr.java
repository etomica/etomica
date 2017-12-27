/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr;

import etomica.atom.AtomType;
import etomica.atom.AtomTypeAgentManager;
import etomica.atom.SpeciesAgentManager;
import etomica.box.BoxAgentManager;
import etomica.box.BoxAgentManager.BoxAgentSource;
import etomica.box.BoxCellManager;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.species.ISpecies;
import etomica.util.Arrays;

public abstract class PotentialMasterNbr extends PotentialMaster {

    private static final AtomTypeAgentManager.AgentSource<PotentialArray> atomTypeAgentSource = new AtomTypeAgentManager.AgentSource<PotentialArray>() {
        @Override
        public PotentialArray makeAgent(AtomType type) {
            return new PotentialArray();
        }

        @Override
        public void releaseAgent(PotentialArray agent, AtomType type) {}
    };

    private static final SpeciesAgentManager.AgentSource<PotentialArray> speciesAgentSource = new SpeciesAgentManager.AgentSource<PotentialArray>() {
        @Override
        public PotentialArray makeAgent(ISpecies type) {
            return new PotentialArray();
        }

        @Override
        public void releaseAgent(PotentialArray agent, ISpecies type) {}
    };

    protected final AtomTypeAgentManager<PotentialArray> rangedAgentManager;
    protected final SpeciesAgentManager<PotentialArray> intraAgentManager;
    protected final Simulation simulation;
    protected IPotential[] allPotentials = new IPotential[0];
    protected BoxAgentManager<? extends BoxCellManager> boxAgentManager;

    public PotentialMasterNbr(Simulation sim, BoxAgentManager<? extends BoxCellManager> boxAgentManager) {
        super();
        simulation = sim;
        this.boxAgentManager = boxAgentManager;
        rangedAgentManager = new AtomTypeAgentManager<>(atomTypeAgentSource, this.simulation);
        intraAgentManager = new SpeciesAgentManager<>(speciesAgentSource, this.simulation);

    }
    
    public PotentialGroup makePotentialGroup(int nBody) {
        return new PotentialGroupNbr(nBody);
    }
    
    public void addPotential(IPotentialMolecular potential, ISpecies[] species) {
        if (!(potential instanceof PotentialGroup)) {
            System.err.println("You gave me a concrete molecule potential and I'm very confused now.  I'll pretend like that's OK but don't hold your breath.");
        }
        super.addPotential(potential, species);
    }

    public void potentialAddedNotify(IPotentialAtomic subPotential, PotentialGroup pGroup) {
        super.potentialAddedNotify(subPotential, pGroup);
        AtomType[] atomTypes = pGroup.getAtomTypes(subPotential);
        if (atomTypes == null) {
        	//change
            if (pGroup.nBody() < 2 && subPotential.getRange() == Double.POSITIVE_INFINITY) {
                boolean found = false;
                for (int i=0; i<allPotentials.length; i++) {
                    if (allPotentials[i] == pGroup) {
                        found = true;
                    }
                }
                if (!found) {
                    allPotentials = Arrays.addObject(allPotentials, pGroup);
                }
                //pGroup is PotentialGroupNbr
                //ADDED S
                if(pGroup.nBody() == 1){
                    ISpecies[] parentType = getSpecies(pGroup);
                    intraAgentManager.getAgent(parentType[0]).addPotential(pGroup);
                }
            }
            else {
                //FIXME what to do with this case?  Fail!
                System.err.println("You have a child-potential of a 2-body PotentialGroup or range-dependent potential, but it's not type-based.  Enjoy crashing or fix bug 85");
            }
            return;
        }
        for (int i=0; i<atomTypes.length; i++) {
            addRangedPotential(subPotential,atomTypes[i]);
        }
        addRangedPotentialForTypes(subPotential, atomTypes);
    }

    protected abstract void addRangedPotentialForTypes(IPotentialAtomic subPotential, AtomType[] atomTypes);

    protected void addRangedPotential(IPotentialAtomic potential, AtomType atomType) {

        PotentialArray potentialAtomType = rangedAgentManager.getAgent(atomType);
        potentialAtomType.addPotential(potential);
        boolean found = false;
        for (int i=0; i<allPotentials.length; i++) {
            if (allPotentials[i] == potential) {
                found = true;
            }
        }
        if (!found) {
            allPotentials = Arrays.addObject(allPotentials, potential);
        }
    }
    
    public void removePotential(IPotentialAtomic potential) {
        super.removePotential(potential);
        if (potential.getRange() < Double.POSITIVE_INFINITY) {
            for (PotentialArray potentialArray : this.rangedAgentManager.getAgents().values()) {
                potentialArray.removePotential(potential);
            }
        }
        else if (potential instanceof PotentialGroup) {
            for (PotentialArray potentialArray : this.intraAgentManager.getAgents().values()) {
                potentialArray.removePotential(potential);
            }
        }
        allPotentials = Arrays.removeObject(allPotentials,potential);
    }

    public PotentialArray getRangedPotentials(AtomType atomType) {
        return rangedAgentManager.getAgent(atomType);
    }

    public PotentialArray getIntraPotentials(ISpecies atomType) {
        return intraAgentManager.getAgent(atomType);
    }

    public final BoxAgentManager<? extends BoxCellManager> getCellAgentManager() {
        return boxAgentManager;
    }
}
