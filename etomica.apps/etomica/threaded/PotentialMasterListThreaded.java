package etomica.threaded;

import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomGroup;
import etomica.atom.AtomPositionDefinition;
import etomica.atom.iterator.IteratorDirective;
import etomica.nbr.PotentialGroupNbr;
import etomica.nbr.cell.PhaseAgentSourceCellManager;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.phase.Phase;
import etomica.phase.PhaseAgentManager;
import etomica.potential.Potential;
import etomica.potential.PotentialArray;
import etomica.potential.PotentialCalculation;
import etomica.space.Space;
import etomica.util.Debug;

public class PotentialMasterListThreaded extends PotentialMasterList {

	PotentialMasterListWorker[] threads;
	int ready1;
	int ready2;
	
	
	public PotentialMasterListThreaded(Space space) {
		super(space);
		// TODO Auto-generated constructor stub
	}

	public PotentialMasterListThreaded(Space space, double range) {
		super(space, range);
		// TODO Auto-generated constructor stub
	}

	public PotentialMasterListThreaded(Space space, double range,
			AtomPositionDefinition positionDefinition) {
		super(space, range, positionDefinition);
		// TODO Auto-generated constructor stub
	}

	public PotentialMasterListThreaded(Space space, double range,
			PhaseAgentSourceCellManager phaseAgentSource) {
		super(space, range, phaseAgentSource);
		// TODO Auto-generated constructor stub
	}

	public PotentialMasterListThreaded(Space space, double range,
			PhaseAgentSourceCellManager phaseAgentSource,
			PhaseAgentManager agentManager) {
		super(space, range, phaseAgentSource, agentManager);
		// TODO Auto-generated constructor stub
	}
	
    public void calculate(Phase phase, IteratorDirective id, PotentialCalculation pc) {
        if(!enabled) return;
        Atom targetAtom = id.getTargetAtom();
        NeighborListManager neighborManager = (NeighborListManager)neighborListAgentManager.getAgent(phase);

        if (targetAtom == null) {
            //no target atoms specified -- do one-target algorithm to SpeciesMaster
            if (Debug.ON && id.direction() != IteratorDirective.Direction.UP) {
                throw new IllegalArgumentException("When there is no target, iterator directive must be up");
            }
            // invoke setPhase on all potentials
            for (int i=0; i<allPotentials.length; i++) {
                allPotentials[i].setPhase(phase);
            }
            
            if(pc instanceof PotentialCalculationThreaded){
            	calculateThreaded(phase, id, (PotentialCalculationThreaded)pc, neighborManager);
            }
            else{
            	//method of super class
            	super.calculate(phase, id, pc);
            }
            
        }
        else {
            //first walk up the tree looking for 1-body range-independent potentials that apply to parents
            Atom parentAtom = targetAtom.getParentGroup();
            while (parentAtom.getType().getDepth() > 1) {
                PotentialArray potentialArray = getIntraPotentials(parentAtom.getType());
                Potential[] potentials = potentialArray.getPotentials();
                for(int i=0; i<potentials.length; i++) {
                    potentials[i].setPhase(phase);
                    ((PotentialGroupNbr)potentials[i]).calculateRangeIndependent(parentAtom,id,pc);
                }
                parentAtom = parentAtom.getParentGroup();
            }                
            PotentialArray potentialArray = (PotentialArray)rangedAgentManager.getAgent(targetAtom.getType());
            Potential[] potentials = potentialArray.getPotentials();
            for(int i=0; i<potentials.length; i++) {
                potentials[i].setPhase(phase);
            }
            calculate(targetAtom, id, pc, neighborManager);
        }
       
        if(lrcMaster != null) {
            lrcMaster.calculate(phase, id, pc);
        }
    }

    protected void calculateThreaded(Phase phase, IteratorDirective id, PotentialCalculationThreaded pc, NeighborListManager neighborManager) {

        //cannot use AtomIterator field because of recursive call
        AtomArrayList list = phase.getSpeciesMaster().getAgentList();
        int size = list.size();
        for (int i=0; i<size; i++) {
            Atom a = list.get(i);
            calculateThreaded(a, id, pc, neighborManager);//recursive call
        }
        pc.writeData();
    }
        
    protected void calculateThreaded(Atom atom, IteratorDirective id, PotentialCalculationThreaded pc, NeighborListManager neighborManager) {
           
        AtomArrayList list = ((AtomGroup)atom).getChildList();
        int size = list.size();
        
    	
			                            
            for(int i=0; i<threads.length; i++){
                synchronized(threads[i]){
                threads[i].startAtom = (i*size)/threads.length;
                threads[i].stopAtom = ((i+1)*size)/threads.length;
                threads[i].threadList = list;
                threads[i].id = id;
                threads[i].pc = pc.getPotentialCalculations()[i];
                threads[i].neighborManager = neighborManager;
                threads[i].greenLight = true;
                threads[i].finished = false;
                threads[i].notifyAll();
               
                }
            }
	
            
			// All threads are running
            
		
		
        // Waiting for all threads to complete, threads report "Finished!"
		for(int i=0; i<threads.length; i++){
			synchronized(threads[i]){
				try{
                   
					if (!threads[i].finished){
                      
                        threads[i].wait();	
					}
				}
				catch(InterruptedException e){
				}
		
			}
		}
		
        
    }
	
	public void setNumThreads(int t){
		threads = new PotentialMasterListWorker[t];
		
		for (int i=0; i<t; i++){
			threads[i] = new PotentialMasterListWorker(i, rangedAgentManager, this);
			threads[i].start();
		}
	}
}
