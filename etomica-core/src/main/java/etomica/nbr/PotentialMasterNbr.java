/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr;

import etomica.atom.AtomType;
import etomica.atom.AtomTypeAgentManager;
import etomica.species.SpeciesAgentManager;
import etomica.box.BoxAgentManager;
import etomica.box.BoxCellManager;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.species.ISpecies;
import etomica.util.Arrays;

import java.util.ArrayList;
import java.util.List;

public abstract class PotentialMasterNbr extends PotentialMaster {

    private final PotentialArray[] rangedPotentials;
    private final PotentialArray[] intraPotentials;
    protected final Simulation simulation;
    protected IPotential[] allPotentials = new IPotential[0];
    protected BoxAgentManager<? extends BoxCellManager> boxAgentManager;

    public PotentialMasterNbr(Simulation sim, BoxAgentManager<? extends BoxCellManager> boxAgentManager) {
        super();
        simulation = sim;
        this.boxAgentManager = boxAgentManager;
        rangedPotentials = sim.getSpeciesList().stream()
                .flatMap(species -> species.getAtomTypes().stream())
                .map(atomType -> new PotentialArray())
                .toArray(PotentialArray[]::new);

        intraPotentials = new PotentialArray[sim.getSpeciesList().size()];
        for (int i = 0; i < intraPotentials.length; i++) {
            intraPotentials[i] = new PotentialArray();
        }
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
                    intraPotentials[parentType[0].getIndex()].addPotential(pGroup);
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

        PotentialArray potentialAtomType = rangedPotentials[atomType.getIndex()];
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
            for (PotentialArray potentialArray : this.rangedPotentials) {
                potentialArray.removePotential(potential);
            }
        }
        else if (potential instanceof PotentialGroup) {
            for (PotentialArray potentialArray : this.intraPotentials) {
                potentialArray.removePotential(potential);
            }
        }
        allPotentials = Arrays.removeObject(allPotentials,potential);
    }

    public final PotentialArray[] getRangedPotentials() {
        return rangedPotentials;
    }

    public final PotentialArray getRangedPotentials(AtomType atomType) {
        return rangedPotentials[atomType.getIndex()];
    }

    public final PotentialArray getIntraPotentials(ISpecies species) {
        return intraPotentials[species.getIndex()];
    }

    public final BoxAgentManager<? extends BoxCellManager> getCellAgentManager() {
        return boxAgentManager;
    }
}
