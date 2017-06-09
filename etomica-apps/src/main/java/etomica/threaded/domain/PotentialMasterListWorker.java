/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.threaded.domain;

import etomica.api.IPotential;
import etomica.api.IPotentialAtomic;
import etomica.atom.*;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.nbr.list.NeighborListManager;
import etomica.potential.IteratorDirective;
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
                IAtom a = threadList.getAtom(i);
                calculate(a, id, pc);
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
            
            IPotentialAtomic potentialThread = ((PotentialThreaded)potentials[i]).getPotentials()[threadNumber];
            
            
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
                    IAtomList list = neighborManager.getUpList(atom)[i];
                    atomPair.atom0 = atom;
                    for (int j=0; j<list.getAtomCount(); j++) {
                        atomPair.atom1 = list.getAtom(j);
                        pc.doCalculation(atomPair, potentialThread);
                    }
                }
                
                // This won't work efficiently with threads
                if (direction != IteratorDirective.Direction.UP) {
                    IAtomList list = neighborManager.getDownList(atom)[i];
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
                    IAtomList list = neighborManager.getUpList(atom)[i];
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
	}
	
	protected void doNBodyStuff(IAtom atom, PotentialCalculation pc, int potentialIndex, IPotentialAtomic potential) {
	        AtomArrayList arrayList = atomsetArrayList.getArrayList();
	        arrayList.clear();
	        arrayList.add(atom);
	        IAtomList[] list = neighborManager.getUpList(atom);
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
    
    public IAtomList threadList;
    protected final int threadNumber;
    public IteratorDirective id;
    public PotentialCalculation pc;
    public NeighborListManager neighborManager;
    public AtomArrayList[][] neighborLists;
    
	//	 things needed for N-body potentials
	protected AtomListWrapper atomsetArrayList;
		
	
}

