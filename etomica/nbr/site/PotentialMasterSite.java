package etomica.nbr.site;

import etomica.atom.AtomAddressManager;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomGroup;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.atom.iterator.AtomsetIteratorPDT;
import etomica.atom.iterator.AtomsetIteratorSinglet;
import etomica.atom.iterator.IteratorDirective;
import etomica.nbr.CriterionAll;
import etomica.nbr.CriterionInterMolecular;
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
import etomica.simulation.Simulation;
import etomica.util.Arrays;
import etomica.util.Debug;

public class PotentialMasterSite extends PotentialMasterNbr {

	/**
	 * Invokes superclass constructor, specifying IteratorFactoryCell
     * for generating molecule iterators.  Sets default nCells of 10 and
     * position definition to null, so that atom type's definition is used
     * to assign cells. 
	 */
	public PotentialMasterSite(Simulation sim, int nCells) {
        this(sim, new PhaseAgentSiteManager(nCells));
    }
    
    public PotentialMasterSite(Simulation sim, PhaseAgentSource phaseAgentSource) {
        this(sim, phaseAgentSource, new PhaseAgentManager(phaseAgentSource));
    }
    
    public PotentialMasterSite(Simulation sim, PhaseAgentSource phaseAgentSource, PhaseAgentManager agentManager) {
        this(sim, phaseAgentSource, agentManager, new Api1ASite(sim.getSpace().D(),agentManager));
    }
    
    protected PotentialMasterSite(Simulation sim, PhaseAgentSource phaseAgentSource, 
            PhaseAgentManager agentManager, AtomsetIteratorPDT neighborIterator) {
        super(sim, phaseAgentSource, agentManager);
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
    public void setCellRange(int newCellRange) {
        cellRange = newCellRange;
    }
    
    protected void addRangedPotentialForTypes(Potential potential, AtomType[] atomType) {
        NeighborCriterion criterion;
        if (atomType.length == 2) {
            criterion = new CriterionTypePair(new CriterionAll(), atomType[0], atomType[1]);
            if (atomType[0].getDepth() > AtomAddressManager.MOLECULE_DEPTH &&
                    atomType[1].getDepth() > AtomAddressManager.MOLECULE_DEPTH) {
                AtomType moleculeType0 = atomType[0].getParentType();
                while (moleculeType0.getDepth() > AtomAddressManager.MOLECULE_DEPTH) {
                    moleculeType0 = moleculeType0.getParentType();
                }
                AtomType moleculeType1 = atomType[1].getParentType();
                while (moleculeType1.getDepth() > AtomAddressManager.MOLECULE_DEPTH) {
                    moleculeType1 = moleculeType1.getParentType();
                }
                if (moleculeType0 == moleculeType1) {
                    criterion = new CriterionInterMolecular(criterion);
                }
            }
        }
        else if (atomType.length == 1) {
            criterion = new CriterionType(new CriterionAll(), atomType[0]);
        }
        else {
            criterion = new CriterionTypesMulti(new CriterionAll(), atomType);
        }
        for (int i=0; i<atomType.length; i++) {
            ((PotentialArray)rangedAgentManager.getAgent(atomType[i])).setCriterion(potential, criterion);
        }
        criteriaArray = (NeighborCriterion[]) Arrays.addObject(criteriaArray, criterion);
    }
    
    /**
     * Returns the criterion used by to determine what atoms interact with the
     * given potential.
     */
    public NeighborCriterion getCriterion(Potential potential) {
        rangedPotentialIterator.reset();
        while (rangedPotentialIterator.hasNext()) {
            PotentialArray potentialArray = (PotentialArray)rangedPotentialIterator.next();
            Potential[] potentials = potentialArray.getPotentials();
            for (int j=0; j<potentials.length; j++) {
                if (potentials[j] == potential) {
                    return potentialArray.getCriteria()[j];
                }
            }
        }
        return null;
    }
    
    /**
     * Sets the criterion associated with the given potential, overriding the 
     * default provided by the PotentialMasterCell.  The criterion can be 
     * configured by calling getCriterion(Potential) and changing the 
     * criterion.  The potential passed to this method must be a potential 
     * handled by this instance.
     */
    public void setCriterion(Potential potential, NeighborCriterion criterion) {
        rangedPotentialIterator.reset();
        while (rangedPotentialIterator.hasNext()) {
            PotentialArray potentialArray = (PotentialArray)rangedPotentialIterator.next();
            Potential[] potentials = potentialArray.getPotentials();
            for (int j=0; j<potentials.length; j++) {
                if (potentials[j] == potential) {
                    potentialArray.setCriterion(potential, criterion);
                    return;
                }
            }
        }
        throw new IllegalArgumentException("Potential "+potential+" is not associated with this PotentialMasterList");
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
        IAtom targetAtom = id.getTargetAtom();
        neighborIterator.setPhase(phase);
        if (targetAtom == null) {
            if (Debug.ON && id.direction() != IteratorDirective.Direction.UP) {
                throw new IllegalArgumentException("When there is no target, iterator directive must be up");
            }
            neighborIterator.setDirection(IteratorDirective.Direction.UP);
            // invoke setPhase on all potentials
            for (int i=0; i<allPotentials.length; i++) {
                allPotentials[i].setPhase(phase);
            }

            //no target atoms specified
            //call calculate with each SpeciesAgent
            AtomSet list = phase.getSpeciesMaster().getAgentList();
            int size = list.getAtomCount();
            for (int i=0; i<size; i++) {
                IAtom a = list.getAtom(i);
                calculate(a, id, pc);//call calculate with the SpeciesAgent
            }
        }
        else {
            // one target atom
            neighborIterator.setDirection(id.direction());
            //first walk up the tree looking for 1-body range-independent potentials that apply to parents
            IAtom pseudoTargetAtom = targetAtom;
            while (pseudoTargetAtom.getType().getDepth() > AtomAddressManager.SPECIES_DEPTH) {
                pseudoTargetAtom = pseudoTargetAtom.getParentGroup();
                PotentialArray potentialArray = getIntraPotentials(pseudoTargetAtom.getType());
                Potential[] potentials = potentialArray.getPotentials();
                for(int i=0; i<potentials.length; i++) {
                    potentials[i].setPhase(phase);
                    ((PotentialGroupNbr)potentials[i]).calculateRangeIndependent(pseudoTargetAtom,id,pc);
                }
            }
            PotentialArray potentialArray = (PotentialArray)rangedAgentManager.getAgent(targetAtom.getType());
            Potential[] potentials = potentialArray.getPotentials();
            for(int i=0; i<potentials.length; i++) {
                potentials[i].setPhase(phase);
            }
            calculate(targetAtom, id, pc);
        }
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
	protected void calculate(IAtom atom, IteratorDirective id, PotentialCalculation pc) {
        PotentialArray potentialArray = (PotentialArray)rangedAgentManager.getAgent(atom.getType());
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
                for (AtomPair pair = (AtomPair)neighborIterator.next(); pair != null;
                     pair = (AtomPair)neighborIterator.next()) {
                    if (nbrCriterion.accept(pair)) {
                        singletPairIterator.setAtom(pair);
                        pc.doCalculation(singletPairIterator, id, p2);
                    }
                }
                break;
            }
        }
            
		//if atom has children, repeat process with them
		if(atom instanceof IAtomGroup) {
            potentialArray = getIntraPotentials(atom.getType());
            potentials = potentialArray.getPotentials();
            for(int i=0; i<potentials.length; i++) {
                ((PotentialGroupNbr)potentials[i]).calculateRangeIndependent(atom,id,pc);
            }

            //cannot use AtomIterator field because of recursive call
            AtomSet list = ((IAtomGroup)atom).getChildList();
            int size = list.getAtomCount();
            for (int i=0; i<size; i++) {
                IAtom a = list.getAtom(i);
                calculate(a, id, pc);//recursive call
            }
		}
	}
    
    private static final long serialVersionUID = 1L;
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
