package etomica.nbr.list;

import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomPositionDefinition;
import etomica.atom.AtomSet;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.SpeciesRoot;
import etomica.atom.iterator.ApiInnerFixed;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.atom.iterator.AtomsetIterator;
import etomica.atom.iterator.AtomsetIteratorBasisDependent;
import etomica.atom.iterator.AtomsetIteratorDirectable;
import etomica.atom.iterator.IteratorDirective;
import etomica.nbr.PotentialMasterNbr;
import etomica.nbr.PotentialCalculationUpdateTypeList.PotentialAtomTypeWrapper;
import etomica.nbr.cell.NeighborCellManager;
import etomica.nbr.cell.PhaseAgentSourceCellManager;
import etomica.phase.Phase;
import etomica.phase.PhaseAgentManager;
import etomica.potential.Potential;
import etomica.potential.PotentialArray;
import etomica.potential.PotentialCalculation;
import etomica.space.Space;

/**
 * PotentialMaster used to implement neighbor listing.  Instance of this
 * class is given as an argument to the Simulation constructor.
 * Criteria specifying whether two atoms are neighbors for a particular potential
 * are specified in the setSpecies method of this class.
 * <br>
 */
public class PotentialMasterList extends PotentialMasterNbr {

	/**
     * Default constructor uses range of 1.0.
     */
    public PotentialMasterList(Space space) {
        this(space,1.0);
    }
    
    /**
     * Constructor specifying space and range for neighbor listing; uses null AtomPositionDefinition.
     */
    public PotentialMasterList(Space space, double range) {
        this(space, range, (AtomPositionDefinition)null);
    }
    
    /**
     * Constructs class using given position definition for all atom cell
     * assignments.
     * 
     * @param positionDefinition
     *            if null, specifies use of atom type's position definition
     */
    public PotentialMasterList(Space space, double range, AtomPositionDefinition positionDefinition) {
        this(space, range, new PhaseAgentSourceCellManager(positionDefinition));
    }
    
    public PotentialMasterList(Space space, double range, PhaseAgentSourceCellManager phaseAgentSource) {
        this(space, range, phaseAgentSource, new PhaseAgentManager(phaseAgentSource,null));
    }
    
    public PotentialMasterList(Space space, double range, PhaseAgentSourceCellManager phaseAgentSource, PhaseAgentManager agentManager) {
        super(space, phaseAgentSource, agentManager);
        neighborManager = new NeighborListManager(this, range, agentManager);
        atomIterator = new AtomIteratorArrayListSimple();
        singletIterator = new AtomIteratorSinglet();
        pairIterator = new ApiInnerFixed(singletIterator, atomIterator);
        swappedPairIterator = new ApiInnerFixed(singletIterator, atomIterator, true);
        cellRange = 2;
    }
    
    public void setRange(double range) {
        neighborManager.setRange(range);
        ((PhaseAgentSourceCellManager)phaseAgentSource).setRange(range);
    }
    
    public double getRange() {
        return neighborManager.getRange();
    }
    
    protected void addToPotentialTypeList(PotentialAtomTypeWrapper wrapper) {
        neighborManager.addCriterion(wrapper.potential.getCriterion(),wrapper.atomTypes);
        super.addToPotentialTypeList(wrapper);
    }

    /**
     * Overrides superclass method to enable direct neighbor-list iteration
     * instead of iteration via species/potential hierarchy. If no target atoms are
     * specified in directive, neighborlist iteration is begun with
     * speciesMaster of phase, and repeated recursively down species hierarchy;
     * if one atom is specified, neighborlist iteration is performed on it and
     * down species hierarchy from it; if two or more atoms are specified,
     * superclass method is invoked.
     */
    public void calculate(Phase phase, IteratorDirective id, PotentialCalculation pc) {
        if(!enabled) return;
        AtomSet targetAtoms = id.getTargetAtoms();
        neighborManager.setPhase(phase);
        if (targetAtoms.count() == 0) {
    		//no target atoms specified -- do one-target algorithm to SpeciesMaster
            PotentialArray potentialArray = getPotentials(phase.getSpeciesMaster().type);
    		    calculate(phase.getSpeciesMaster(), idUp, pc, potentialArray.getPotentials(), potentialArray.getIterators());
    		    if(lrcMaster != null) {
    		        lrcMaster.calculate(phase, id, pc);
    		    }
    	    }
        else if (targetAtoms instanceof Atom) {
    		// one target atom
            PotentialArray potentialArray = getPotentials(((Atom)targetAtoms).type);
            calculate((Atom)targetAtoms, id, pc, potentialArray.getPotentials(), potentialArray.getIterators());
            if(lrcMaster != null) {
                lrcMaster.calculate(phase, id, pc);
            }
        }
    	    else {
    	        //more than one target atom
    	        super.calculate(phase, id, pc);
    	    }
    }//end calculate
	
    /**
     * Performs given PotentialCalculation using potentials/neighbors associated
     * with the given atom (if any).  Then, if atom is not a leaf atom, iteration over
     * child atoms is performed and process is repeated (recursively) with each on down
     * the hierarchy until leaf atoms are reached.
     */
    //TODO make a "TerminalGroup" node that permits child atoms but indicates that no potentials apply directly to them
    private void calculate(Atom atom, IteratorDirective id, PotentialCalculation pc, Potential[] potentials, AtomsetIterator[] iterators) {
        singletIterator.setAtom(atom);
        IteratorDirective.Direction direction = id.direction();
        for(int i=0; i<potentials.length; i++) {
            if (!potentials[i].getCriterion().isRangeDependent()) {
                // not range-dependent, so assume intragroup!  let's hope we're right.
                
                // if you're hitting a ClassCastExcpetion here, it's because you didn't 
                // give your potential an appropriate criterion
                ((AtomsetIteratorBasisDependent)iterators[i]).setBasis(atom.node.parentGroup());
                ((AtomsetIteratorBasisDependent)iterators[i]).setTarget(atom);
                ((AtomsetIteratorDirectable)iterators[i]).setDirection(direction);
                pc.doCalculation(iterators[i],id,potentials[i]);
                continue;
            }
            switch (potentials[i].nBody()) {
            case 1:
                boolean[] potential1BodyArray = neighborManager.getPotential1BodyList(atom).getInteractingList();
                if (i < potential1BodyArray.length && potential1BodyArray[i]) {
                    pc.doCalculation(singletIterator, id, potentials[i]);
                }
                break;
            case 2:
                AtomArrayList[] list;
                if (direction != IteratorDirective.Direction.DOWN) {
                    list = neighborManager.getUpList(atom);
//                  list.length may be less than potentials.length, if atom hasn't yet interacted with another using one of the potentials
                    if(i < list.length) {
                        atomIterator.setList(list[i]);
                        //System.out.println("Up :"+atomIterator.size());
                        pc.doCalculation(pairIterator, id, potentials[i]);
                    }
                }
                if (direction != IteratorDirective.Direction.UP) {
                    list = neighborManager.getDownList(atom);
                    if(i < list.length) {
                        atomIterator.setList(list[i]);
                        //System.out.println("Dn :"+atomIterator.size());
                        pc.doCalculation(swappedPairIterator, id, potentials[i]);
                    }
                }
                break;//switch
            }//end of switch
        }//end of for
        
        //if atom has children, repeat process with them
        if(!atom.node.isLeaf()) {
            //cannot use AtomIterator field because of recursive call
            AtomArrayList list = ((AtomTreeNodeGroup) atom.node).childList;
            int size = list.size();
            for (int i=0; i<size; i++) {
                Atom a = list.get(i);
                PotentialArray potentialArray = getPotentials(a.type);
                Potential[] childPotentials = potentialArray.getPotentials();
                AtomsetIterator[] childIterators = potentialArray.getIterators();
                calculate(a, id, pc, childPotentials, childIterators);//recursive call
            }
        }
    }

    public NeighborListManager getNeighborManager() {return neighborManager;}

    
    public NeighborCellManager getNbrCellManager(Phase phase) {
        phaseAgentManager.setRoot((SpeciesRoot)phase.getSpeciesMaster().node.parentGroup());
        NeighborCellManager[] cellManagers = (NeighborCellManager[])phaseAgentManager.getAgents();
        NeighborCellManager manager = cellManagers[phase.getIndex()];
        manager.setPotentialRange(neighborManager.getRange());
        manager.setCellRange(cellRange);
        return manager;
    }

    public void setCellRange(int newCellRange) {
        cellRange = newCellRange;
    }

    public int getCellRange() {
        return cellRange;
    }
    
    private final AtomIteratorArrayListSimple atomIterator;
    private final AtomIteratorSinglet singletIterator;
    private final ApiInnerFixed pairIterator;
    private final ApiInnerFixed swappedPairIterator;
    private final NeighborListManager neighborManager;
    private int cellRange;
    private final IteratorDirective idUp = new IteratorDirective();
    
}
