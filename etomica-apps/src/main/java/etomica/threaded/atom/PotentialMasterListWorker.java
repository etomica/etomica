/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.threaded.atom;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.IPotential;
import etomica.api.IPotentialAtomic;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomListWrapper;
import etomica.atom.AtomPair;
import etomica.atom.AtomSetSinglet;
import etomica.atom.AtomTypeAgentManager;
import etomica.atom.MoleculeArrayList;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.atom.iterator.IteratorDirective;
import etomica.nbr.PotentialGroupNbr;
import etomica.nbr.list.NeighborListManager;
import etomica.potential.PotentialArray;
import etomica.potential.PotentialCalculation;
import etomica.threaded.PotentialThreaded;

/*
 * XXX: This class is unlikely to work properly due to changes made since this
 * class was last used.
 */
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
            for(int i=startAtom; i<stopAtom; i++){
                IMolecule molecule = threadList.getMolecule(i-startAtom);
                PotentialArray potentialArray = pmlt.getIntraPotentials(molecule.getType());
                IPotential[] potentials = potentialArray.getPotentials();
                for(int j=0; j<potentials.length; j++) {
                
                    // Extracts thread-specific potential for intra-molecular atoms
                    IPotential potentialIntraThread = ((PotentialThreaded)potentials[j]).getPotentials()[threadNumber];
                    ((PotentialGroupNbr)potentialIntraThread).calculateRangeIndependent(molecule,id.direction(),null,pc);
                }
                
                //cannot use AtomIterator field because of recursive call
                IAtomList list = molecule.getChildList();
                int size = list.getAtomCount();
                for (int j=0; j<size; j++) {
                    IAtom a = list.getAtom(j);
                    calculate(a, id, pc);//recursive call
                }
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
                    IAtomList list = neighborLists[atom.getIndex()-startAtom][i];
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
                    IAtomList list = neighborLists[atom.getIndex()+startAtom][i];
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
    
    public void fillNeighborListArray(int threadNumber, int numThreads, NeighborListManager nm, Box box){
        
        // Make reference to neighbor lists
        IMoleculeList list = box.getMoleculeList();
        int size = list.getMoleculeCount();
        
        int startAtom = (threadNumber*size)/numThreads;
        int stopAtom = ((threadNumber+1)*size)/numThreads;
        
        threadList = new MoleculeArrayList();
        
        this.neighborManager = nm;
        neighborManager.updateLists();
        neighborLists = new IAtomList[stopAtom-startAtom][];
                
        for(int i=0; i<(stopAtom-startAtom); i++){
            threadList.add(list.getMolecule(i+startAtom));
            neighborLists[i] = neighborManager.getUpList(list.getMolecule(i+startAtom).getChildList().getAtom(0));
            
        }
            
    
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
    
    public MoleculeArrayList threadList;
    public int startAtom;
    public int stopAtom;
    protected final int threadNumber;
    public IteratorDirective id;
    public PotentialCalculation pc;
    public NeighborListManager neighborManager;
    public IAtomList[][] neighborLists;
    
	//	 things needed for N-body potentials
	protected AtomListWrapper atomsetArrayList;
		
	
}

