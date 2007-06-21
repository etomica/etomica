package etomica.threaded.domain;

import etomica.atom.AtomArrayList;
import etomica.atom.AtomSet;
import etomica.atom.AtomTypeAgentManager;
import etomica.atom.AtomsetArrayList;
import etomica.atom.IAtom;
import etomica.atom.IAtomGroup;
import etomica.atom.iterator.ApiInnerFixed;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.atom.iterator.AtomsetIteratorSinglet;
import etomica.atom.iterator.IteratorDirective;
import etomica.nbr.PotentialGroupNbr;
import etomica.nbr.list.NeighborListManager;
import etomica.potential.IPotential;
import etomica.potential.PotentialArray;
import etomica.potential.PotentialCalculation;
import etomica.threaded.PotentialThreaded;

public class PotentialMasterListWorker extends Thread {

	public PotentialMasterListWorker (int threadNumber, AtomTypeAgentManager rangedAgentManager, PotentialMasterListThreaded pmlt){
		
        
        
        this.pmlt = pmlt;
        this.threadNumber = threadNumber;
		this.rangedAgentManager = rangedAgentManager;
		atomIterator = new AtomIteratorArrayListSimple();
        singletIterator = new AtomIteratorSinglet();
        pairIterator = new ApiInnerFixed(singletIterator, atomIterator);
        swappedPairIterator = new ApiInnerFixed(singletIterator, atomIterator, true);
		
		
	}

	public void run(){
		
		while (active){
			
			// FIRST Wait on MAIN thread
			synchronized(this){
				try{
					if 	(!active){ break;}
                    
                    // Reset finished flag
                    
                   
                    // Thread is waiting for instructions from PotentialMasterListThreaded
					if (!greenLight){
						
                        
                        wait();	
					}
                  
				}
				catch(InterruptedException e){
					if(!active){ break;}
				}
              
			}
			
            
            
                        
           threadCalculate = System.currentTimeMillis(); 
            // Thread completes objective
            for(int i=0; i<threadList.getAtomCount(); i++){
                calculate(threadList.getAtom(i), id, pc);       
            }
           threadCalculate = System.currentTimeMillis()-threadCalculate;
           
        
        
            System.out.println("Thread"+threadNumber+" start/stop is "+threadCalculate);          
                      
           
                      
			synchronized(this){
			    // Stops thread before work is started again
                greenLight = false;
                
                finished = true;
				notifyAll();		
			}
					
		}
		
		
	}
	
	
	public void calculate(IAtom atom, IteratorDirective id, PotentialCalculation pc) {
		
		 singletIterator.setAtom(atom);
	     IteratorDirective.Direction direction = id.direction();
	     PotentialArray potentialArray = (PotentialArray)rangedAgentManager.getAgent(atom.getType());
	     IPotential[] potentials = potentialArray.getPotentials();
		
		for(int i=0; i<potentials.length; i++) {
            
            IPotential potentialThread = ((PotentialThreaded)potentials[i]).getPotentials()[threadNumber];
            
            
            switch (potentialThread.nBody()) {
            case 1:
                boolean[] potential1BodyArray = neighborManager.getPotential1BodyList(atom).getInteractingList();
                if (potential1BodyArray[i]) {
                    pc.doCalculation(singletIterator, id, potentialThread);
                }
                break;
            case 2:
                AtomSet[] list;
                if (direction != IteratorDirective.Direction.DOWN) {
                    list = neighborManager.getUpList(atom);
//                  list.length may be less than potentials.length, if atom hasn't yet interacted with another using one of the potentials
                    atomIterator.setList(list[i]);
                    //System.out.println("Up :"+atomIterator.size());
                    
                   // threadCalculate2 -= System.nanoTime();
                    pc.doCalculation(pairIterator, id, potentialThread);
                   // threadCalculate2 += System.nanoTime();
                }
                
                // This won't work efficiently with threads
                if (direction != IteratorDirective.Direction.UP) {
                    list = neighborManager.getDownList(atom);
                    atomIterator.setList(list[i]);
                    //System.out.println("Dn :"+atomIterator.size());
                    pc.doCalculation(swappedPairIterator, id, potentialThread);
                }
                break;//switch
            case Integer.MAX_VALUE: //N-body
                // instantiate lazily so other simulations don't have to carry this stuff around
                if (atomsetArrayList == null) {
                    atomsetArrayList = new AtomsetArrayList();
                    singletSetIterator = new AtomsetIteratorSinglet(atomsetArrayList);
                }
                // do the calculation considering the current Atom as the 
                // "central" Atom.
                doNBodyStuff(atom, id, pc, i, potentialThread);
                if (direction != IteratorDirective.Direction.UP) {
                    // must have a target and be doing "both"
                    // we have to do the calculation considering each of the 
                    // target's neighbors
                    list = neighborManager.getUpList(atom);
                    if (i < list.length) {
                        AtomSet iList = list[i];
                        for (int j=0; j<iList.getAtomCount(); j++) {
                            IAtom otherAtom = iList.getAtom(j);
                            doNBodyStuff(otherAtom, id, pc, i, potentialThread);
                        }
                    }
                    list = neighborManager.getDownList(atom);
                    if (i < list.length) {
                        AtomSet iList = list[i];
                        for (int j=0; j<iList.getAtomCount(); j++) {
                            IAtom otherAtom = iList.getAtom(j);
                            doNBodyStuff(otherAtom, id, pc, i, potentialThread);
                        }
                    }
                }
                
            }//end of switch
        }//end of for

		//		if atom has children, repeat process with them
        if(atom instanceof IAtomGroup) {
            potentialArray = pmlt.getIntraPotentials(atom.getType());
            potentials = potentialArray.getPotentials();
            for(int i=0; i<potentials.length; i++) {
            
                // Extracts thread-specific potential for intra-molecular atoms
                IPotential potentialIntraThread = ((PotentialThreaded)potentials[i]).getPotentials()[threadNumber];
                ((PotentialGroupNbr)potentialIntraThread).calculateRangeIndependent(atom,id,pc);
            }
            
            //cannot use AtomIterator field because of recursive call
            AtomSet list = ((IAtomGroup)atom).getChildList();
            int size = list.getAtomCount();
            for (int i=0; i<size; i++) {
                IAtom a = list.getAtom(i);
                calculate(a, id, pc);//recursive call
            }
        }
	}
	
	protected void doNBodyStuff(IAtom atom, IteratorDirective id, PotentialCalculation pc, int potentialIndex, IPotential potential) {
	        AtomArrayList arrayList = atomsetArrayList.getArrayList();
	        arrayList.clear();
	        arrayList.add(atom);
	        AtomSet[] list = neighborManager.getUpList(atom);
	        if (potentialIndex < list.length) {
	            arrayList.addAll(list[potentialIndex]);
	        }
	        list = neighborManager.getDownList(atom);
	        if (potentialIndex < list.length) {
	            arrayList.addAll(list[potentialIndex]);
	        }
	        pc.doCalculation(singletSetIterator, id, potential);
	        arrayList.clear();
    }
    
        
    protected final PotentialMasterListThreaded pmlt;
	protected final AtomIteratorArrayListSimple atomIterator;
	protected final AtomIteratorSinglet singletIterator;
	protected final ApiInnerFixed pairIterator;
	protected final ApiInnerFixed swappedPairIterator;
	protected final AtomTypeAgentManager rangedAgentManager;
    
    protected boolean active = true;
    protected boolean greenLight = false;
    protected boolean finished = false;
    
    public double threadCalculate;
   // public double threadCalculate2;
    
    public AtomArrayList threadList;
    protected final int threadNumber;
    public IteratorDirective id;
    public PotentialCalculation pc;
    public NeighborListManager neighborManager;
    public AtomArrayList[][] neighborLists;
    
	//	 things needed for N-body potentials
	protected AtomsetArrayList atomsetArrayList;
	protected AtomsetIteratorSinglet singletSetIterator;
		
	
}

