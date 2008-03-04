package etomica.threaded.domain;

import etomica.atom.AtomPositionDefinition;
import etomica.atom.AtomTypeMolecule;
import etomica.api.IAtom;
import etomica.api.IBox;
import etomica.atom.IAtomLeaf;
import etomica.api.IMolecule;
import etomica.atom.iterator.IteratorDirective;
import etomica.box.BoxAgentManager;
import etomica.lattice.CellLattice;
import etomica.nbr.PotentialGroupNbr;
import etomica.nbr.cell.BoxAgentSourceCellManager;
import etomica.nbr.cell.Cell;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.api.IPotential;
import etomica.potential.PotentialArray;
import etomica.potential.PotentialCalculation;
import etomica.api.ISimulation;
import etomica.space.Space;
import etomica.threaded.IPotentialCalculationThreaded;
import etomica.util.Debug;

public class PotentialMasterListThreaded extends PotentialMasterList {

    private static final long serialVersionUID = 1L;
    PotentialMasterListWorker[] threads;
	BoxAgentManager agentManagerThreaded;
	
	
	public PotentialMasterListThreaded(ISimulation sim, Space _space) {
		this(sim, 0, _space);
		// TODO Auto-generated constructor stub
	}

	public PotentialMasterListThreaded(ISimulation sim, double range, Space _space) {
        this(sim, range, (AtomPositionDefinition)null, _space);
		// TODO Auto-generated constructor stub
	}

	public PotentialMasterListThreaded(ISimulation sim, double range,
			AtomPositionDefinition positionDefinition, Space _space) {
        this(sim, range, new BoxAgentSourceCellManager(positionDefinition, _space), _space);
		// TODO Auto-generated constructor stub
	}

	public PotentialMasterListThreaded(ISimulation sim, double range,
			BoxAgentSourceCellManager boxAgentSource, Space _space) {
		this(sim, range, boxAgentSource, new BoxAgentManager(boxAgentSource), _space);
		// TODO Auto-generated constructor stub
	}

	public PotentialMasterListThreaded(ISimulation sim, double range,
			BoxAgentSourceCellManager boxAgentSource,
			BoxAgentManager agentManager, Space _space) {
		super(sim, range, boxAgentSource, agentManager, new NeighborListAgentSourceThreaded(range, _space));
        agentManagerThreaded = new BoxAgentManager(new BoxAgentSourceCellManagerThreaded(null, _space), sim, true);
	}
	
    public NeighborCellManagerThreaded getNbrCellManagerThreaded(IBox box) {
        return (NeighborCellManagerThreaded)agentManagerThreaded.getAgent(box);
    }
    
    public void calculate(IBox box, IteratorDirective id, PotentialCalculation pc) {
        if(!enabled) return;
        IAtom targetAtom = id.getTargetAtom();
        NeighborListManager neighborManager = (NeighborListManager)neighborListAgentManager.getAgent(box);

        if (targetAtom == null) {
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
            if (targetAtom instanceof IAtomLeaf) {
                IMolecule molecule = ((IAtomLeaf)targetAtom).getParentGroup();
                PotentialArray potentialArray = getIntraPotentials((AtomTypeMolecule)molecule.getType());
                IPotential[] potentials = potentialArray.getPotentials();
                for(int i=0; i<potentials.length; i++) {
                    potentials[i].setBox(box);
                    ((PotentialGroupNbr)potentials[i]).calculateRangeIndependent(molecule,id,pc);
                }
                potentialArray = (PotentialArray)rangedAgentManager.getAgent(targetAtom.getType());
                potentials = potentialArray.getPotentials();
                for(int i=0; i<potentials.length; i++) {
                    potentials[i].setBox(box);
                }
                calculate((IAtomLeaf)targetAtom, id, pc, neighborManager);
            }
            else {
                PotentialArray potentialArray = (PotentialArray)rangedAgentManager.getAgent(targetAtom.getType());
                IPotential[] potentials = potentialArray.getPotentials();
                for(int i=0; i<potentials.length; i++) {
                    potentials[i].setBox(box);
                }
                calculate((IMolecule)targetAtom, id, pc, neighborManager);
            }
        }
       
        if(lrcMaster != null) {
            lrcMaster.calculate(box, id, pc);
        }
    }

    protected void calculateThreaded(IBox box, IteratorDirective id, IPotentialCalculationThreaded pc, NeighborListManager neighborManager) {
    	
			                            
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
	
	public void setNumThreads(int t, IBox box){
        
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
            this.space = _space;
        }

        public Class getAgentClass() {
            return NeighborListManagerThreaded.class;
        }

        public Object makeAgent(IBox box) {
            return new NeighborListManagerThreaded((PotentialMasterListThreaded)potentialMaster, range, box, space);
        }

        private final Space space;
        private static final long serialVersionUID = 1L;

      
    }



}
