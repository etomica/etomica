package etomica.threaded.domain;

import etomica.api.IAtom;
import etomica.api.IAtomSet;
import etomica.api.IMolecule;
import etomica.api.IPotential;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomPair;
import etomica.atom.AtomSetSinglet;
import etomica.atom.AtomTypeAgentManager;
import etomica.atom.AtomTypeMolecule;
import etomica.atom.AtomsetArrayList;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.atom.iterator.IteratorDirective;
import etomica.nbr.PotentialGroupNbr;
import etomica.nbr.list.NeighborListManager;
import etomica.potential.PotentialArray;
import etomica.potential.PotentialCalculation;
import etomica.threaded.PotentialThreaded;

public class PotentialMasterListWorker extends Thread {

	public PotentialMasterListWorker (int threadNumber, AtomTypeAgentManager rangedAgentManager, PotentialMasterListThreaded pmlt){
		
        
        
        this.pmlt = pmlt;
        this.threadNumber = threadNumber;
		this.rangedAgentManager = rangedAgentManager;
		atomIterator = new AtomIteratorArrayListSimple();
		atomSetSinglet = new AtomSetSinglet();
		atomPair = new AtomPair();
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
		
	     IteratorDirective.Direction direction = id.direction();
	     PotentialArray potentialArray = (PotentialArray)rangedAgentManager.getAgent(atom.getType());
	     IPotential[] potentials = potentialArray.getPotentials();
		
		for(int i=0; i<potentials.length; i++) {
            
            IPotential potentialThread = ((PotentialThreaded)potentials[i]).getPotentials()[threadNumber];
            
            
            switch (potentialThread.nBody()) {
            case 1:
                boolean[] potential1BodyArray = neighborManager.getPotential1BodyList(atom).getInteractingList();
                if (potential1BodyArray[i]) {
                    atomSetSinglet.atom = atom;
                    pc.doCalculation(atomSetSinglet, potentialThread);
                }
                break;
            case 2:
                if (direction != IteratorDirective.Direction.DOWN) {
                    IAtomSet list = neighborManager.getUpList(atom)[i];
                    atomPair.atom0 = atom;
                    for (int j=0; j<list.getAtomCount(); j++) {
                        atomPair.atom1 = list.getAtom(j);
                        pc.doCalculation(atomPair, potentialThread);
                    }
                }
                
                // This won't work efficiently with threads
                if (direction != IteratorDirective.Direction.UP) {
                    IAtomSet list = neighborManager.getDownList(atom)[i];
                    atomPair.atom1 = atom;
                    for (int j=0; j<list.getAtomCount(); j++) {
                        atomPair.atom0 = list.getAtom(j);
                        pc.doCalculation(atomPair, potentialThread);
                    }
                }
                break;//switch
            case Integer.MAX_VALUE: //N-body
                // do the calculation considering the current Atom as the 
                // "central" Atom.
                doNBodyStuff(atom, pc, i, potentialThread);
                if (direction != IteratorDirective.Direction.UP) {
                    // must have a target and be doing "both"
                    // we have to do the calculation considering each of the 
                    // target's neighbors
                    IAtomSet list = neighborManager.getUpList(atom)[i];
                    for (int j=0; j<list.getAtomCount(); j++) {
                        IAtom otherAtom = list.getAtom(j);
                        doNBodyStuff(otherAtom, pc, i, potentialThread);
                    }
                    list = neighborManager.getDownList(atom)[i];
                    for (int j=0; j<list.getAtomCount(); j++) {
                        IAtom otherAtom = list.getAtom(j);
                        doNBodyStuff(otherAtom, pc, i, potentialThread);
                    }
                }
                
            }//end of switch
        }//end of for

		//		if atom has children, repeat process with them
        if(atom instanceof IMolecule) {
            potentialArray = pmlt.getIntraPotentials((AtomTypeMolecule)atom.getType());
            potentials = potentialArray.getPotentials();
            for(int i=0; i<potentials.length; i++) {
            
                // Extracts thread-specific potential for intra-molecular atoms
                IPotential potentialIntraThread = ((PotentialThreaded)potentials[i]).getPotentials()[threadNumber];
                ((PotentialGroupNbr)potentialIntraThread).calculateRangeIndependent((IMolecule)atom,id,pc);
            }
            
            //cannot use AtomIterator field because of recursive call
            IAtomSet list = ((IMolecule)atom).getChildList();
            int size = list.getAtomCount();
            for (int i=0; i<size; i++) {
                IAtom a = list.getAtom(i);
                calculate(a, id, pc);//recursive call
            }
        }
	}
	
	protected void doNBodyStuff(IAtom atom, PotentialCalculation pc, int potentialIndex, IPotential potential) {
	        AtomArrayList arrayList = atomsetArrayList.getArrayList();
	        arrayList.clear();
	        arrayList.add(atom);
	        IAtomSet[] list = neighborManager.getUpList(atom);
	        if (potentialIndex < list.length) {
	            arrayList.addAll(list[potentialIndex]);
	        }
	        list = neighborManager.getDownList(atom);
	        if (potentialIndex < list.length) {
	            arrayList.addAll(list[potentialIndex]);
	        }
	        pc.doCalculation(atomsetArrayList, potential);
	        arrayList.clear();
    }
    
        
    protected final PotentialMasterListThreaded pmlt;
	protected final AtomIteratorArrayListSimple atomIterator;
	protected final AtomSetSinglet atomSetSinglet;
	protected final AtomPair atomPair;
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
		
	
}

