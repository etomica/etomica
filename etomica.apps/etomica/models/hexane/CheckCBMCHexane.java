package etomica.models.hexane;

import etomica.action.Action;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomGroup;
import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.phase.Phase;
import etomica.space.IVector;

/**
 * 
 * @author cribbin
 *
 */
public class CheckCBMCHexane implements Action {

    public CheckCBMCHexane(Phase p){
        phase = p;
        energyMeter = new MeterPotentialEnergy(phase.getSimulation().getPotentialMaster());
        energyMeter.setPhase(phase);
        
        moleculeIterator = new AtomIteratorAllMolecules(phase);
        moleculeIterator.reset();

        atomList = new AtomArrayList(6);
        
        vex = phase.getSpace().makeVector();
        booink = -1;
    }
    
    public void actionPerformed() {
        if(energyMeter.getDataAsScalar() != 0.0){
            throw new RuntimeException("Bad Zoot!  Non-zero potential energy!");
        }
    
        moleculeIterator.reset();
        
        while (moleculeIterator.hasNext()){
            atomList = ((AtomGroup)moleculeIterator.nextAtom()).getChildList();
            double length;
            double tol = 0.000005;
            
            for(int i = 0 ; i < atomList.size() - 1; i++){
                vex.E(((AtomLeaf)atomList.get(i)).getCoord().getPosition());
                vex.ME(((AtomLeaf)atomList.get(i+1)).getCoord().getPosition());
                length = Math.sqrt(vex.squared());
                length -= 0.4;
                if(length > tol){
                    throw new RuntimeException("Bad Zoot!  The bond lengh changed!");
                }
            }
            
        }
     
        System.out.println(booink++);
    }

    MeterPotentialEnergy energyMeter;
    Phase phase;
    AtomIteratorAllMolecules moleculeIterator;
    AtomArrayList atomList;
    IVector vex;
    int booink;
}