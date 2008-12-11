package etomica.dimer;

import etomica.api.IAtom;
import etomica.api.IBox;
import etomica.api.IMoleculeList;
import etomica.api.ISimulation;
import etomica.api.ISpecies;
import etomica.atom.iterator.IteratorDirective;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.PotentialCalculation;
import etomica.space.ISpace;
import etomica.util.Debug;

public class PotentialMasterListDimer extends PotentialMasterList{

    public PotentialMasterListDimer(ISimulation sim, ISpace space) {
        super(sim, space);
        
    }
    
    public void setSpecies(ISpecies[] species){
        this.species = species;
    }
    
   public void calculate(IBox box, IteratorDirective id, PotentialCalculation pc) {
        if(!enabled) return;
        IAtom targetAtom = id.getTargetAtom();
        if (targetAtom != null) {
            super.calculate(box, id, pc);
            return;
        }
        NeighborListManager neighborManager = (NeighborListManager)neighborListAgentManager.getAgent(box);
        
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
