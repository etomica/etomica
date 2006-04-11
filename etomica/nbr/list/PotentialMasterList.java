package etomica.nbr.list;

import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomPositionDefinition;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.AtomType;
import etomica.atom.AtomsetArrayList;
import etomica.atom.SpeciesRoot;
import etomica.atom.iterator.ApiInnerFixed;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.atom.iterator.AtomsetIteratorSinglet;
import etomica.atom.iterator.IteratorDirective;
import etomica.nbr.PotentialGroupNbr;
import etomica.nbr.PotentialMasterNbr;
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
    
    protected void addRangedPotentialToList(Potential potential, AtomType atomType) {
        neighborManager.addCriterion(potential.getCriterion(),atomType);
        super.addRangedPotentialToList(potential, atomType);
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
        Atom targetAtom = id.getTargetAtom();
        neighborManager.setPhase(phase);
        if (targetAtom == null) {
            //no target atoms specified -- do one-target algorithm to SpeciesMaster
            targetAtom = phase.getSpeciesMaster();
            id = idUp;
            // invoke setPhase on all potentials
            for (int i=0; i<allPotentials.length; i++) {
                allPotentials[i].setPhase(phase);
            }
        }
        else {
            //first walk up the tree looking for 1-body range-independent potentials that apply to parents
            Atom parentAtom = targetAtom.node.parentGroup();
            while (parentAtom.type.getDepth() > 2) {
                PotentialArray potentialArray = getIntraPotentials(parentAtom.type);
                Potential[] potentials = potentialArray.getPotentials();
                for(int i=0; i<potentials.length; i++) {
                    potentials[i].setPhase(phase);
                    ((PotentialGroupNbr)potentials[i]).calculateRangeIndependent(parentAtom,id,pc);
                }
                parentAtom = parentAtom.node.parentGroup();
            }                
            PotentialArray potentialArray = getRangedPotentials(targetAtom.type);
            Potential[] potentials = potentialArray.getPotentials();
            for(int i=0; i<potentials.length; i++) {
                potentials[i].setPhase(phase);
            }
        }
        calculate(targetAtom, id, pc);
        if(lrcMaster != null) {
            lrcMaster.calculate(phase, id, pc);
        }
    }
	
    /**
     * Performs given PotentialCalculation using potentials/neighbors associated
     * with the given atom (if any).  Then, if atom is not a leaf atom, iteration over
     * child atoms is performed and process is repeated (recursively) with each on down
     * the hierarchy until leaf atoms are reached.
     */
    //TODO make a "TerminalGroup" node that permits child atoms but indicates that no potentials apply directly to them
    private void calculate(Atom atom, IteratorDirective id, PotentialCalculation pc) {
        singletIterator.setAtom(atom);
        IteratorDirective.Direction direction = id.direction();
        PotentialArray potentialArray = getRangedPotentials(atom.type);
        Potential[] potentials = potentialArray.getPotentials();
        for(int i=0; i<potentials.length; i++) {
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
            case Integer.MAX_VALUE: //N-body
                // instantiate lazily so other simulations don't have to carry this stuff around
                if (atomsetArrayList == null) {
                    atomsetArrayList = new AtomsetArrayList();
                    singletSetIterator = new AtomsetIteratorSinglet(atomsetArrayList);
                }
                // do the calculation considering the current Atom as the 
                // "central" Atom.
                doNBodyStuff(atom, id, pc, i, potentials[i]);
                if (direction != IteratorDirective.Direction.UP) {
                    // must have a target and be doing "both"
                    // we have to do the calculation considering each of the 
                    // target's neighbors
                    list = neighborManager.getUpList(atom);
                    if (i < list.length) {
                        AtomArrayList iList = list[i];
                        for (int j=0; j<iList.size(); j++) {
                            Atom otherAtom = iList.get(j);
                            doNBodyStuff(otherAtom, id, pc, i, potentials[i]);
                        }
                    }
                    list = neighborManager.getDownList(atom);
                    if (i < list.length) {
                        AtomArrayList iList = list[i];
                        for (int j=0; j<iList.size(); j++) {
                            Atom otherAtom = iList.get(j);
                            doNBodyStuff(otherAtom, id, pc, i, potentials[i]);
                        }
                    }
                }
                
            }//end of switch
        }//end of for
        
        //if atom has children, repeat process with them
        if(!atom.node.isLeaf()) {
            potentialArray = getIntraPotentials(atom.type);
            potentials = potentialArray.getPotentials();
            for(int i=0; i<potentials.length; i++) {
                ((PotentialGroupNbr)potentials[i]).calculateRangeIndependent(atom,id,pc);
            }
            
            //cannot use AtomIterator field because of recursive call
            AtomArrayList list = ((AtomTreeNodeGroup) atom.node).childList;
            int size = list.size();
            for (int i=0; i<size; i++) {
                Atom a = list.get(i);
                calculate(a, id, pc);//recursive call
            }
        }
    }
    
    /**
     * Invokes the PotentialCalculation for the given Atom with its up and down
     * neighbors as a single AtomSet.
     */
    protected void doNBodyStuff(Atom atom, IteratorDirective id, PotentialCalculation pc, int potentialIndex, Potential potential) {
        AtomArrayList arrayList = atomsetArrayList.getArrayList();
        arrayList.clear();
        arrayList.add(atom);
        AtomArrayList[] list = neighborManager.getUpList(atom);
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
    
    // things needed for N-body potentials
    private AtomsetArrayList atomsetArrayList;
    private AtomsetIteratorSinglet singletSetIterator;
}
