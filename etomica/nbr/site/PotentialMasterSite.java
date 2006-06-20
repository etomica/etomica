/*
 * History
 * Created on Sep 20, 2004 by kofke
 */
package etomica.nbr.site;

import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomPair;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.AtomType;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.atom.iterator.AtomsetIteratorPDT;
import etomica.atom.iterator.AtomsetIteratorSinglet;
import etomica.atom.iterator.IteratorDirective;
import etomica.nbr.CriterionAll;
import etomica.nbr.CriterionSimple;
import etomica.nbr.CriterionType;
import etomica.nbr.CriterionTypePair;
import etomica.nbr.CriterionTypesMulti;
import etomica.nbr.NeighborCriterion;
import etomica.nbr.PotentialGroupNbr;
import etomica.nbr.PotentialMasterNbr;
import etomica.phase.Phase;
import etomica.phase.PhaseAgentManager;
import etomica.phase.PhaseAgentManager.PhaseAgentSource;
import etomica.potential.Potential;
import etomica.potential.Potential2;
import etomica.potential.PotentialArray;
import etomica.potential.PotentialCalculation;
import etomica.space.Space;
import etomica.util.Arrays;
import etomica.util.Debug;

public class PotentialMasterSite extends PotentialMasterNbr {

	/**
	 * Invokes superclass constructor, specifying IteratorFactoryCell
     * for generating molecule iterators.  Sets default nCells of 10 and
     * position definition to null, so that atom type's definition is used
     * to assign cells. 
	 */
	public PotentialMasterSite(Space space, int nCells) {
        this(space, new PhaseAgentSiteManager(nCells));
    }
    
    public PotentialMasterSite(Space space, PhaseAgentSource phaseAgentSource) {
        this(space, phaseAgentSource, new PhaseAgentManager(phaseAgentSource,null));
    }
    
    public PotentialMasterSite(Space space, PhaseAgentSource phaseAgentSource, PhaseAgentManager agentManager) {
        this(space, phaseAgentSource, agentManager, new Api1ASite(space.D(),agentManager));
    }
    
    protected PotentialMasterSite(Space space, PhaseAgentSource phaseAgentSource, 
            PhaseAgentManager agentManager, AtomsetIteratorPDT neighborIterator) {
        super(space, phaseAgentSource, agentManager);
        singletAtomIterator = new AtomIteratorSinglet();
		singletPairIterator = new AtomsetIteratorSinglet(2);
        this.neighborIterator = neighborIterator;
	}
    
    /**
     * @return Returns the cellRange.
     */
    public int getCellRange() {
        return cellRange;
    }

    /**
     * @param cellRange The cellRange to set.
     */
    public void setCellRange(int cellRange) {
        this.cellRange = cellRange;
    }
    
    protected void addRangedPotentialForTypes(Potential potential, AtomType[] atomType) {
        NeighborCriterion criterion;
        if (atomType.length == 2) {
            criterion = new CriterionTypePair(new CriterionAll(), atomType[0], atomType[1]);
        }
        else if (atomType.length == 1) {
            criterion = new CriterionType(new CriterionAll(), atomType[0]);
        }
        else {
            criterion = new CriterionTypesMulti(new CriterionAll(), atomType);
        }
        for (int i=0; i<atomType.length; i++) {
            rangedPotentialAtomTypeList[atomType[i].getIndex()].setCriterion(potential, criterion);
        }
        criteriaArray = (NeighborCriterion[]) Arrays.addObject(criteriaArray, criterion);
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
        if (!enabled)
            return;
        for (int i=0; i<criteriaArray.length; i++) {
            criteriaArray[i].setPhase(phase);
        }
        Atom targetAtom = id.getTargetAtom();
        neighborIterator.setPhase(phase);
        if (targetAtom == null) {
            //no target atoms specified -- do one-target algorithm to SpeciesMaster
            targetAtom = phase.getSpeciesMaster();
            if (Debug.ON && id.direction() != IteratorDirective.Direction.UP) {
                throw new IllegalArgumentException("When there is no target, iterator directive must be up");
            }
            neighborIterator.setDirection(IteratorDirective.Direction.UP);
            // invoke setPhase on all potentials
            for (int i=0; i<allPotentials.length; i++) {
                allPotentials[i].setPhase(phase);
            }
        }
        else {
            // one target atom
            neighborIterator.setDirection(id.direction());
            //first walk up the tree looking for 1-body range-independent potentials that apply to parents
            Atom pseudoTargetAtom = targetAtom;
            while (pseudoTargetAtom.type.getDepth() > 3) {
                pseudoTargetAtom = pseudoTargetAtom.node.parentGroup();
                PotentialArray potentialArray = getIntraPotentials(pseudoTargetAtom.type);
                Potential[] potentials = potentialArray.getPotentials();
                for(int i=0; i<potentials.length; i++) {
                    potentials[i].setPhase(phase);
                    ((PotentialGroupNbr)potentials[i]).calculateRangeIndependent(pseudoTargetAtom,id,pc);
                }
            }                
            PotentialArray potentialArray = getRangedPotentials(targetAtom.type);
            Potential[] potentials = potentialArray.getPotentials();
            for(int i=0; i<potentials.length; i++) {
                potentials[i].setPhase(phase);
            }
        }
        calculate(targetAtom, id, pc);
        if (lrcMaster != null) {
            lrcMaster.calculate(phase, id, pc);
        }
    }
	
    /**
     * Performs given PotentialCalculation using potentials/neighbors associated
     * with the given atom (if any).  Then, if atom is not a leaf atom, iteration over
     * child atoms is performed and process is repeated (recursively) with each on down
     * the hierarchy until leaf atoms are reached.
     */
    //TODO make a "TerminalGroup" indicator in type that permits child atoms but indicates that no potentials apply directly to them
	protected void calculate(Atom atom, IteratorDirective id, PotentialCalculation pc) {
        PotentialArray potentialArray = getRangedPotentials(atom.type);
        Potential[] potentials = potentialArray.getPotentials();
        NeighborCriterion[] criteria = potentialArray.getCriteria();

        for(int i=0; i<potentials.length; i++) {
            switch (potentials[i].nBody()) {
            case 1:
                singletAtomIterator.setAtom(atom);
                pc.doCalculation(singletAtomIterator, id, potentials[i]);
                break;
            case 2:
                Potential2 p2 = (Potential2) potentials[i];
                NeighborCriterion nbrCriterion = criteria[i];
                neighborIterator.setTarget(atom);
                neighborIterator.reset();
                while (neighborIterator.hasNext()) {
                    AtomPair pair = (AtomPair)neighborIterator.next();
                    if (nbrCriterion.accept(pair)) {
                        singletPairIterator.setAtom(pair);
                        pc.doCalculation(singletPairIterator, id, p2);
                    }
                }
                break;
            }
        }
            
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
    
    private final AtomIteratorSinglet singletAtomIterator;
	private final AtomsetIteratorSinglet singletPairIterator;
    private int cellRange;
    protected final AtomsetIteratorPDT neighborIterator;
    private NeighborCriterion[] criteriaArray = new NeighborCriterion[0];
    
    public static class PhaseAgentSiteManager implements PhaseAgentSource {
        public PhaseAgentSiteManager(int nCells) {
            this.nCells = nCells;
        }
        
        public Class getAgentClass() {
            return NeighborSiteManager.class;
        }
        
        public Object makeAgent(Phase phase) {
            return new NeighborSiteManager(phase,nCells);
        }
        
        public void releaseAgent(Object agent) {
        }
        
        private final int nCells;
    }
}
