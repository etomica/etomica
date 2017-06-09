/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.threaded.atom;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.atom.IMolecule;
import etomica.atom.IMoleculeList;
import etomica.api.IPotential;
import etomica.simulation.Simulation;
import etomica.atom.IMoleculePositionDefinition;
import etomica.atom.iterator.IteratorDirective;
import etomica.box.BoxAgentManager;
import etomica.nbr.PotentialGroupNbr;
import etomica.nbr.list.BoxAgentSourceCellManagerList;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.PotentialArray;
import etomica.potential.PotentialCalculation;
import etomica.space.Space;
import etomica.threaded.IPotentialCalculationThreaded;
import etomica.util.Debug;

public class PotentialMasterListThreaded extends PotentialMasterList {

    private static final long serialVersionUID = 1L;
    PotentialMasterListWorker[] threads;
	int ready1;
	int ready2;
	
	
	public PotentialMasterListThreaded(Simulation sim, Space _space) {
		super(sim, _space);
	}

	public PotentialMasterListThreaded(Simulation sim, Space _space, double range) {
		super(sim, range, _space);
	}

	public PotentialMasterListThreaded(Simulation sim, Space _space, double range,
                                       IMoleculePositionDefinition positionDefinition) {
		super(sim, range, positionDefinition, _space);
	}

	public PotentialMasterListThreaded(Simulation sim, Space _space, double range,
                                       BoxAgentSourceCellManagerList boxAgentSource) {
		super(sim, range, boxAgentSource, _space);
	}

	public PotentialMasterListThreaded(Simulation sim, Space _space, double range,
                                       BoxAgentSourceCellManagerList boxAgentSource,
                                       BoxAgentManager agentManager) {
		super(sim, range, boxAgentSource, agentManager, _space);
	}
	
    public void calculate(Box box, IteratorDirective id, PotentialCalculation pc) {
        if(!enabled) return;
        IAtom targetAtom = id.getTargetAtom();
        IMolecule targetMolecule = id.getTargetMolecule();
        NeighborListManager neighborManager = neighborListAgentManager.getAgent(box);

        if (targetAtom == null && targetMolecule == null) {
            //no target atoms specified -- do one-target algorithm to SpeciesMaster
            if (Debug.ON && id.direction() != IteratorDirective.Direction.UP) {
                throw new IllegalArgumentException("When there is no target, iterator directive must be up");
            }
            // invoke setBox on all potentials
            for (int i=0; i<allPotentials.length; i++) {
                allPotentials[i].setBox(box);
            }
            
            if(pc instanceof IPotentialCalculationThreaded){
            	calculateThreaded(box, id, (IPotentialCalculationThreaded)pc, neighborManager);
            }
            else{
            	//method of super class
            	super.calculate(box, id, pc);
            }
            
        }
        else {
            //first walk up the tree looking for 1-body range-independent potentials that apply to parents
            if (targetAtom != null) {
                IMolecule molecule = targetAtom.getParentGroup();
                PotentialArray potentialArray = getIntraPotentials(molecule.getType());
                IPotential[] potentials = potentialArray.getPotentials();
                for(int i=0; i<potentials.length; i++) {
                    potentials[i].setBox(box);
                    ((PotentialGroupNbr)potentials[i]).calculateRangeIndependent(molecule,id.direction(),targetAtom,pc);
                }
                
                potentialArray = (PotentialArray)rangedAgentManager.getAgent(targetAtom.getType());
                potentials = potentialArray.getPotentials();
                for(int i=0; i<potentials.length; i++) {
                    potentials[i].setBox(box);
                }
                calculate(targetAtom, id.direction(), pc, neighborManager);
            }
            else {
                PotentialArray potentialArray = (PotentialArray)rangedAgentManager.getAgent(targetAtom.getType());
                IPotential[] potentials = potentialArray.getPotentials();
                for(int i=0; i<potentials.length; i++) {
                    potentials[i].setBox(box);
                }
                calculate(targetMolecule, id.direction(), pc, neighborManager);
            }
        }
       
        if(lrcMaster != null) {
            lrcMaster.calculate(box, id, pc);
        }
    }

    protected void calculateThreaded(Box box, IteratorDirective id, IPotentialCalculationThreaded pc, NeighborListManager neighborManager) {

        //cannot use AtomIterator field because of recursive call
        IMoleculeList list = box.getMoleculeList();
        int size = list.getMoleculeCount();
			                            
            for(int i=0; i<threads.length; i++){
                synchronized(threads[i]){
                threads[i].startAtom = (i*size)/threads.length;
                threads[i].stopAtom = ((i+1)*size)/threads.length;
                threads[i].id = id;
                threads[i].pc = pc.getPotentialCalculations()[i];
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
		
        pc.writeData();
        
    }
	
	public void setNumThreads(int t, Box box){
		threads = new PotentialMasterListWorker[t];
				
        for (int i=0; i<t; i++){
			threads[i] = new PotentialMasterListWorker(i, rangedAgentManager, this);
			threads[i].fillNeighborListArray(i, t, neighborListAgentManager.getAgent(box), box);
            threads[i].start();
		}
           
	}
}
