package etomica.nbr.list;

import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomPositionDefinition;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.AtomType;
import etomica.atom.AtomsetArrayList;
import etomica.atom.iterator.ApiInnerFixed;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.atom.iterator.AtomsetIteratorSinglet;
import etomica.atom.iterator.IteratorDirective;
import etomica.nbr.CriterionAdapter;
import etomica.nbr.CriterionInterMolecular;
import etomica.nbr.CriterionSimple;
import etomica.nbr.CriterionType;
import etomica.nbr.CriterionTypePair;
import etomica.nbr.CriterionTypesMulti;
import etomica.nbr.NeighborCriterion;
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
import etomica.util.Debug;

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
    
    /**
     * Sets the range that determines how far to look for neighbors.  This 
     * range must be greater than the range of the longest-range potential.
     */
    public void setRange(double range) {
        if (range <= 0) {
            throw new IllegalArgumentException("Range must be greater than 0");
        }
        neighborManager.setRange(range);
        ((PhaseAgentSourceCellManager)phaseAgentSource).setRange(range);
        recomputeCriteriaRanges();
    }
    
    /**
     * Returns the range that determines how far to look for neighbors.
     */
    public double getRange() {
        return neighborManager.getRange();
    }
    
    /**
     * Sets the safety factor.  The default value, 0.4, is appropriate for most
     * simulations.  Valid values range from 0 to 0.5 (non-inclusive).  The 
     * safety factor determines how far an atom can travel from its original
     * position before the neighbor lists are reconstructed.  0.5 means the 
     * atom can travel half of its neighbor range.  If another atom also 
     * travels half-way then the Atoms could interact without the 
     * PotentialMaster recognizing they are neighbors.
     * 
     * High values of the safetyFactor make it more probable that an Atom will
     * move too far between checks and interact without the PotentialMaster 
     * knowing.  Smaller values make it less probable, but slow down the 
     * simulation due to more frequent neighbor list constructing.
     */
    public void setSafetyFactor(double newSafetyFactor) {
        if (newSafetyFactor <= 0 || newSafetyFactor >= 0.5) {
            throw new IllegalArgumentException("Safety factor must be between 0 and 0.5");
        }
        safetyFactor = newSafetyFactor;
        recomputeCriteriaRanges();
    }
    
    
    /**
     * Returns the safety factor.
     */
    public double getSafetyFactor() {
        return safetyFactor;
    }

    /**
     * Adds the potential as a ranged potential that applies to the given 
     * AtomTypes.  This method creates a criterion for the potential and 
     * notifies the NeighborListManager of its existence.
     */
    protected void addRangedPotentialForTypes(Potential potential, AtomType[] atomType) {
        // we'll fix the neighbor range later in recomputeCriteriaRanges
        // 0 guarantees the simulation to be hosed if our range is less than the potential range
        // (since recomputeCriteriaRange will bail in that case)
        CriterionSimple rangedCriterion = new CriterionSimple(getSimulation(), potential.getRange(), 0.0);
        NeighborCriterion criterion;
        if (atomType.length == 2) {
            criterion = new CriterionTypePair(rangedCriterion, atomType[0], atomType[1]);
            if (atomType[0].getDepth() > 3 && atomType[1].getDepth() > 3) {
                AtomType moleculeType0 = atomType[0].getParentType();
                while (moleculeType0.getDepth() > 3) {
                    moleculeType0 = moleculeType0.getParentType();
                }
                AtomType moleculeType1 = atomType[1].getParentType();
                while (moleculeType1.getDepth() > 3) {
                    moleculeType1 = moleculeType1.getParentType();
                }
                if (moleculeType0 == moleculeType1) {
                    criterion = new CriterionInterMolecular(criterion);
                }
            }
        }
        else if (atomType.length == 1) {
            criterion = new CriterionType(rangedCriterion, atomType[0]);
        }
        else {
            criterion = new CriterionTypesMulti(rangedCriterion, atomType);
        }
        neighborManager.addCriterion(criterion);
        for (int i=0; i<atomType.length; i++) {
            rangedPotentialAtomTypeList[atomType[i].getIndex()].setCriterion(potential, criterion);
        }
        if (potential.getRange() > maxPotentialRange) {
            maxPotentialRange = potential.getRange();
        }
        recomputeCriteriaRanges();
    }
    
    /**
     * Recomputes the range for all criterion based on our own range.  The 
     * range for each criterion is set so that they all have an equal max
     * displacement.  Our nominal neighbor range is used for the criterion
     * with the longest potential range.
     */
    protected void recomputeCriteriaRanges() {
        double maxDisplacement = (getRange() - maxPotentialRange) * safetyFactor;
        if (maxDisplacement < 0) {
            // someone probably added a long ranged potential and hasn't updated the PotentialMaster's range
            // when they do, we'll get called again.
            // if they don't, the simulation will probably crash
            return;
        }
        for (int i=0; i<rangedPotentialAtomTypeList.length; i++) {
            NeighborCriterion[] criteria = rangedPotentialAtomTypeList[i].getCriteria();
            Potential[] potentials = rangedPotentialAtomTypeList[i].getPotentials();
            // this will double (or more) count criteria that apply to multiple atom types, but it won't hurt us
            for (int j=0; j<criteria.length; j++) {
                CriterionSimple rangedCriterion = getRangedCriterion(criteria[j]);
                if (rangedCriterion != null) {
                    double newRange = maxDisplacement/safetyFactor + potentials[j].getRange();
                    rangedCriterion.setNeighborRange(newRange);
                    rangedCriterion.setSafetyFactor(safetyFactor);
                }                    
            }
        }
    }
    
    /**
     * Returns the criterion used by to determine what atoms interact with the
     * given potential.
     */
    public NeighborCriterion getCriterion(Potential potential) {
        for (int i=0; i<rangedPotentialAtomTypeList.length; i++) {
            Potential[] potentials = rangedPotentialAtomTypeList[i].getPotentials();
            for (int j=0; j<potentials.length; j++) {
                if (potentials[j] == potential) {
                    return rangedPotentialAtomTypeList[i].getCriteria()[j];
                }
            }
        }
        return null;
    }
    
    /**
     * Convenience method to return the wrapped range-dependent criterion, if 
     * one exists
     */
    private static CriterionSimple getRangedCriterion(NeighborCriterion criterion) {
        if (criterion instanceof CriterionSimple) {
            return (CriterionSimple)criterion;
        }
        if (criterion instanceof CriterionAdapter) {
            return getRangedCriterion(((CriterionAdapter)criterion).getWrappedCriterion());
        }
        return null;
    }
    
    /**
     * Sets the criterion associated with the given potential, overriding the 
     * default provided by the PotentialMasterList.  The criterion can be 
     * configured by calling getCriterion(Potential) and changing the 
     * criterion.  The potential passed to this method must be a potential 
     * handled by this instance.
     */
    public void setCriterion(Potential potential, NeighborCriterion criterion) {
        NeighborCriterion oldCriterion = getCriterion(potential);
        if (oldCriterion != null) {
            neighborManager.removeCriterion(oldCriterion);
        }
        boolean success = false;
        for (int i=0; i<rangedPotentialAtomTypeList.length; i++) {
            Potential[] potentials = rangedPotentialAtomTypeList[i].getPotentials();
            for (int j=0; j<potentials.length; j++) {
                if (potentials[j] == potential) {
                    success = true;
                    rangedPotentialAtomTypeList[i].setCriterion(potential, criterion);
                    break;
                }
            }
        }
        if (success) {
        	neighborManager.addCriterion(criterion);
        	return;
        }
        throw new IllegalArgumentException("Potential "+potential+" is not associated with this PotentialMasterList");
    }
    
    public void removePotential(Potential potential) {
        for (int i=0; i<rangedPotentialAtomTypeList.length; i++) {
            Potential[] potentials = rangedPotentialAtomTypeList[i].getPotentials();
            for (int j=0; j<potentials.length; j++) {
                if (potentials[j] == potential) {
                    // found it!
                    neighborManager.removeCriterion(rangedPotentialAtomTypeList[i].getCriteria()[j]);
                    break;
                }
            }
        }

        super.removePotential(potential);
        
        if (potential.getRange() == maxPotentialRange) {
            maxPotentialRange = 0;
            for (int i=0; i<allPotentials.length; i++) {
                double pRange = allPotentials[i].getRange();
                if (pRange == Double.POSITIVE_INFINITY) {
                    continue;
                }
                if (pRange > maxPotentialRange) {
                    maxPotentialRange = pRange;
                }
            }
            recomputeCriteriaRanges();
        }
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
            if (Debug.ON && id.direction() != IteratorDirective.Direction.UP) {
                throw new IllegalArgumentException("When there is no target, iterator directive must be up");
            }
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
                if (potential1BodyArray[i]) {
                    pc.doCalculation(singletIterator, id, potentials[i]);
                }
                break;
            case 2:
                AtomArrayList[] list;
                if (direction != IteratorDirective.Direction.DOWN) {
                    list = neighborManager.getUpList(atom);
//                  list.length may be less than potentials.length, if atom hasn't yet interacted with another using one of the potentials
                    atomIterator.setList(list[i]);
                    //System.out.println("Up :"+atomIterator.size());
                    pc.doCalculation(pairIterator, id, potentials[i]);
                }
                if (direction != IteratorDirective.Direction.UP) {
                    list = neighborManager.getDownList(atom);
                    atomIterator.setList(list[i]);
                    //System.out.println("Dn :"+atomIterator.size());
                    pc.doCalculation(swappedPairIterator, id, potentials[i]);
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
        PhaseAgentManager phaseAgentManager = getCellAgentManager();
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
    private double maxPotentialRange = 0;
    private double safetyFactor = 0.4;
    
    // things needed for N-body potentials
    private AtomsetArrayList atomsetArrayList;
    private AtomsetIteratorSinglet singletSetIterator;
}
