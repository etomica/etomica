/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.dimer;

import etomica.box.Box;
import etomica.molecule.IMoleculeList;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculation;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.species.ISpecies;
import etomica.util.Debug;

public class PotentialMasterListDimer extends PotentialMasterList{

    public PotentialMasterListDimer(Simulation sim, Space space) {
        super(sim, space);
        
    }
    
    public void setSpecies(ISpecies[] species){
        this.species = species;
    }
    
   public void calculate(Box box, IteratorDirective id, PotentialCalculation pc) {
        if(!enabled) return;
        if (id.getTargetAtom() != null || id.getTargetMolecule() != null) {
            super.calculate(box, id, pc);
            return;
        }
        NeighborListManager neighborManager = neighborListAgentManager.getAgent(box);
        
        if (Debug.ON && id.direction() != IteratorDirective.Direction.UP) {
            throw new IllegalArgumentException("When there is no target, iterator directive must be up");
        }
        // invoke setBox on all potentials
        for (int i=0; i<allPotentials.length; i++) {
            allPotentials[i].setBox(box);
        }

        //no target atoms specified
        //call calculate with each SpeciesAgent
        for(int j=0; j<species.length; j++){    
            IMoleculeList list = box.getMoleculeList(species[j]);
            int size = list.getMoleculeCount();
            for (int i=0; i<size; i++) {
                calculate(list.getMolecule(i), id.direction(), pc, neighborManager);//call calculate with the SpeciesAgent
            }
        }
   }
    
    
    public ISpecies [] species;
}
