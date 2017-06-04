/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr;

import etomica.atom.IAtomType;
import etomica.api.IPotential;
import etomica.api.IPotentialAtomic;
import etomica.api.IPotentialMolecular;
import etomica.simulation.Simulation;
import etomica.api.ISpecies;
import etomica.atom.AtomTypeAgentManager;
import etomica.atom.SpeciesAgentManager;
import etomica.box.BoxAgentManager;
import etomica.box.BoxAgentManager.BoxAgentSource;
import etomica.box.BoxCellManager;
import etomica.potential.PotentialArray;
import etomica.potential.PotentialGroup;
import etomica.potential.PotentialMaster;
import etomica.util.Arrays;

public abstract class PotentialMasterNbr extends PotentialMaster implements AtomTypeAgentManager.AgentSource, SpeciesAgentManager.AgentSource {

    public PotentialMasterNbr(Simulation sim, BoxAgentSource<? extends BoxCellManager> boxAgentSource,
                              BoxAgentManager<? extends BoxCellManager> boxAgentManager) {
        super();
        simulation = sim;
        this.boxAgentSource = boxAgentSource;
        this.boxAgentManager = boxAgentManager;
        rangedAgentManager = new AtomTypeAgentManager(this);
        intraAgentManager = new SpeciesAgentManager(this);

        rangedAgentManager.init(sim);
        intraAgentManager.init(sim);
        rangedPotentialIterator = rangedAgentManager.makeIterator();
        intraPotentialIterator = intraAgentManager.makeIterator();
        boxAgentManager.setSimulation(sim);
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
        IAtomType[] atomTypes = pGroup.getAtomTypes(subPotential);
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
                    allPotentials = (IPotential[])etomica.util.Arrays.addObject(allPotentials, pGroup);
                }
                //pGroup is PotentialGroupNbr
                //ADDED S
                if(pGroup.nBody() == 1){
                    ISpecies[] parentType = getSpecies(pGroup);
                    ((PotentialArray)intraAgentManager.getAgent(parentType[0])).addPotential(pGroup);                	
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

    protected abstract void addRangedPotentialForTypes(IPotentialAtomic subPotential, IAtomType[] atomTypes);
    
    protected void addRangedPotential(IPotentialAtomic potential, IAtomType atomType) {
        
        PotentialArray potentialAtomType = (PotentialArray)rangedAgentManager.getAgent(atomType);
        potentialAtomType.addPotential(potential);
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
    
    public void removePotential(IPotentialAtomic potential) {
        super.removePotential(potential);
        if (potential.getRange() < Double.POSITIVE_INFINITY) {
            rangedPotentialIterator.reset();
            while (rangedPotentialIterator.hasNext()) {
                ((PotentialArray)rangedPotentialIterator.next()).removePotential(potential);
            }
        }
        else if (potential instanceof PotentialGroup) {
            intraPotentialIterator.reset();
            while (intraPotentialIterator.hasNext()) {
                ((PotentialArray)intraPotentialIterator.next()).removePotential(potential);
            }
        }
        allPotentials = (IPotential[])Arrays.removeObject(allPotentials,potential);
    }
    
    public PotentialArray getRangedPotentials(IAtomType atomType) {
        return (PotentialArray)rangedAgentManager.getAgent(atomType);
    }

    public PotentialArray getIntraPotentials(ISpecies atomType) {
        return (PotentialArray)intraAgentManager.getAgent(atomType);
    }
    
    public final BoxAgentManager<? extends BoxCellManager> getCellAgentManager() {
        return boxAgentManager;
    }
    
    public Class getSpeciesAgentClass() {
        return PotentialArray.class;
    }
    
    public Object makeAgent(IAtomType type) {
        return new PotentialArray();
    }
    
    public void releaseAgent(Object agent, IAtomType type) {
    }

    public Object makeAgent(ISpecies type) {
        return new PotentialArray();
    }
    
    public void releaseAgent(Object agent, ISpecies type) {
    }

    /**
     * Returns the simulation associated with this PotentialMaster
     */
    public Simulation getSimulation() {
        return simulation;
    }

    protected AtomTypeAgentManager.AgentIterator rangedPotentialIterator;
    protected SpeciesAgentManager.AgentIterator intraPotentialIterator;
    protected final AtomTypeAgentManager rangedAgentManager;
    protected final SpeciesAgentManager intraAgentManager;
    protected IPotential[] allPotentials = new IPotential[0];
    protected BoxAgentSource<? extends BoxCellManager> boxAgentSource;
    protected final Simulation simulation;
    protected BoxAgentManager<? extends BoxCellManager> boxAgentManager;
}
