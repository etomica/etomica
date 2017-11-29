/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.threaded.domain;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.box.BoxAgentManager;
import etomica.lattice.CellLattice;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculePositionDefinition;
import etomica.nbr.PotentialGroupNbr;
import etomica.nbr.cell.Cell;
import etomica.nbr.cell.NeighborCellManager;
import etomica.nbr.list.BoxAgentSourceCellManagerList;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.IPotential;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialArray;
import etomica.potential.PotentialCalculation;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.threaded.IPotentialCalculationThreaded;
import etomica.util.Debug;

public class PotentialMasterListThreaded extends PotentialMasterList {

    PotentialMasterListWorker[] threads;
	BoxAgentManager agentManagerThreaded;
	
	
	public PotentialMasterListThreaded(Simulation sim, Space _space) {
		this(sim, 0, _space);
	}

	public PotentialMasterListThreaded(Simulation sim, double range, Space _space) {
        this(sim, range, (IMoleculePositionDefinition)null, _space);
	}

	public PotentialMasterListThreaded(Simulation sim, double range,
                                       IMoleculePositionDefinition positionDefinition, Space _space) {
        this(sim, range, new BoxAgentSourceCellManagerList(sim, positionDefinition, _space), _space);
	}

	public PotentialMasterListThreaded(Simulation sim, double range,
                                       BoxAgentSourceCellManagerList boxAgentSource, Space _space) {
        this(sim, range, boxAgentSource, new BoxAgentManager<NeighborCellManager>(boxAgentSource, NeighborCellManager.class, sim), _space);
    }

	public PotentialMasterListThreaded(Simulation sim, double range,
                                       BoxAgentSourceCellManagerList boxAgentSource,
                                       BoxAgentManager<NeighborCellManager> agentManager, Space _space) {
		super(sim, range, boxAgentSource, agentManager, new NeighborListAgentSourceThreaded(range, _space), _space);
        agentManagerThreaded = new BoxAgentManager<NeighborCellManagerThreaded>(new BoxAgentSourceCellManagerThreaded(sim, null, _space), NeighborCellManagerThreaded.class, sim);
	}
	
    public NeighborCellManagerThreaded getNbrCellManagerThreaded(Box box) {
        return (NeighborCellManagerThreaded)agentManagerThreaded.getAgent(box);
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
    	
			                            
            for(int i=0; i<threads.length; i++){
                synchronized(threads[i]){
                threads[i].neighborManager = neighborManager;
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
        
        //Sets the number of domains to the number of threads
		NeighborCellManagerThreaded neighborCellManagerThreaded = (NeighborCellManagerThreaded)agentManagerThreaded.getAgent(box);
        neighborCellManagerThreaded.setNumCells(t);
        
        CellLattice lattice = neighborCellManagerThreaded.getLattice();
        
        
        threads = new PotentialMasterListWorker[t];
        
        for (int i=0; i<t; i++){
			threads[i] = new PotentialMasterListWorker(i, rangedAgentManager, this);
            
            threads[i].threadList = ((Cell)lattice.site(lattice.latticeIndex(i))).occupants();
            threads[i].start();
		}
           
	}


    protected static class NeighborListAgentSourceThreaded extends NeighborListAgentSource{
        
        
        public NeighborListAgentSourceThreaded(double range, Space _space) {

            super(range, _space);
        }

        public NeighborListManagerThreaded makeAgent(Box box) {
            return new NeighborListManagerThreaded((PotentialMasterListThreaded)potentialMaster, range, box, space);
        }
    }



}
