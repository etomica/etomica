package etomica.nbr.site;

import etomica.api.IAtom;
import etomica.api.IAtomSet;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IPotential;
import etomica.api.ISimulation;
import etomica.atom.AtomPair;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.AtomTypeMolecule;
import etomica.atom.IAtomLeaf;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.atom.iterator.AtomsetIteratorPDT;
import etomica.atom.iterator.AtomsetIteratorSinglet;
import etomica.atom.iterator.IteratorDirective;
import etomica.box.BoxAgentManager;
import etomica.box.BoxAgentManager.BoxAgentSource;
import etomica.nbr.CriterionAll;
import etomica.nbr.CriterionInterMolecular;
import etomica.nbr.CriterionType;
import etomica.nbr.CriterionTypePair;
import etomica.nbr.CriterionTypesMulti;
import etomica.nbr.NeighborCriterion;
import etomica.nbr.PotentialGroupNbr;
import etomica.nbr.PotentialMasterNbr;
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
	public PotentialMasterSite(ISimulation sim, int nCells, Space _space) {
        this(sim, new BoxAgentSiteManager(nCells, _space));
    }
    
    public PotentialMasterSite(ISimulation sim, BoxAgentSource boxAgentSource) {
        this(sim, boxAgentSource, new BoxAgentManager(boxAgentSource));
    }
    
    public PotentialMasterSite(ISimulation sim, BoxAgentSource boxAgentSource, BoxAgentManager agentManager) {
        this(sim, boxAgentSource, agentManager, new Api1ASite(sim.getSpace().D(),agentManager));
    }
    
    protected PotentialMasterSite(ISimulation sim, BoxAgentSource boxAgentSource, 
            BoxAgentManager agentManager, AtomsetIteratorPDT neighborIterator) {
        super(sim, boxAgentSource, agentManager);
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
    
    protected void addRangedPotentialForTypes(IPotential potential, IAtomType[] atomType) {
        NeighborCriterion criterion;
        if (atomType.length == 2) {
            criterion = new CriterionTypePair(new CriterionAll(), atomType[0], atomType[1]);
            if ((atomType[0] instanceof AtomTypeLeaf) &&
                    (atomType[1] instanceof AtomTypeLeaf)) {
                IAtomType moleculeType0 = ((AtomTypeLeaf)atomType[0]).getParentType();
                IAtomType moleculeType1 = ((AtomTypeLeaf)atomType[1]).getParentType();
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
    public NeighborCriterion getCriterion(IPotential potential) {
        rangedPotentialIterator.reset();
        while (rangedPotentialIterator.hasNext()) {
            PotentialArray potentialArray = (PotentialArray)rangedPotentialIterator.next();
            IPotential[] potentials = potentialArray.getPotentials();
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
    public void setCriterion(IPotential potential, NeighborCriterion criterion) {
        rangedPotentialIterator.reset();
        while (rangedPotentialIterator.hasNext()) {
            PotentialArray potentialArray = (PotentialArray)rangedPotentialIterator.next();
            IPotential[] potentials = potentialArray.getPotentials();
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
     * speciesMaster of box, and repeated recursively down species hierarchy;
     * if one atom is specified, neighborlist iteration is performed on it and
     * down species hierarchy from it; if two or more atoms are specified,
     * superclass method is invoked.
     */
    public void calculate(IBox box, IteratorDirective id, PotentialCalculation pc) {
        if (!enabled)
            return;
        for (int i=0; i<criteriaArray.length; i++) {
            criteriaArray[i].setBox(box);
        }
        IAtom targetAtom = id.getTargetAtom();
        neighborIterator.setBox(box);
        if (targetAtom == null) {
            if (Debug.ON && id.direction() != IteratorDirective.Direction.UP) {
                throw new IllegalArgumentException("When there is no target, iterator directive must be up");
            }
            neighborIterator.setDirection(IteratorDirective.Direction.UP);
            // invoke setBox on all potentials
            for (int i=0; i<allPotentials.length; i++) {
                allPotentials[i].setBox(box);
            }

            //no target atoms specified
            //call calculate with each SpeciesAgent
            IAtomSet list = box.getMoleculeList();
            int size = list.getAtomCount();
            for (int i=0; i<size; i++) {
                calculate((IMolecule)list.getAtom(i), id, pc);//call calculate with the SpeciesAgent
            }
        }
        else {
            // one target atom
            neighborIterator.setDirection(id.direction());
            PotentialArray potentialArray = (PotentialArray)rangedAgentManager.getAgent(targetAtom.getType());
            IPotential[] potentials = potentialArray.getPotentials();
            for(int i=0; i<potentials.length; i++) {
                potentials[i].setBox(box);
            }
            if (targetAtom instanceof IAtomLeaf) {
                //walk up the tree looking for 1-body range-independent potentials that apply to parents
                IMolecule parentMolecule = ((IAtomLeaf)targetAtom).getParentGroup();
                potentialArray = getIntraPotentials((AtomTypeMolecule)parentMolecule.getType());
                potentials = potentialArray.getPotentials();
                for(int i=0; i<potentials.length; i++) {
                    potentials[i].setBox(box);
                    ((PotentialGroupNbr)potentials[i]).calculateRangeIndependent(parentMolecule,id,pc);
                }
                calculate((IAtomLeaf)targetAtom, id, pc);
            }
            else {
                calculate((IMolecule)targetAtom, id, pc);
            }
        }
        if (lrcMaster != null) {
            lrcMaster.calculate(box, id, pc);
        }
    }
	
    /**
     * Performs given PotentialCalculation using potentials/neighbors associated
     * with the given atom (if any).  Then, if atom is not a leaf atom, iteration over
     * child atoms is performed and process is repeated (recursively) with each on down
     * the hierarchy until leaf atoms are reached.
     */
    //TODO make a "TerminalGroup" indicator in type that permits child atoms but indicates that no potentials apply directly to them
	protected void calculate(IMolecule atom, IteratorDirective id, PotentialCalculation pc) {
        PotentialArray potentialArray = (PotentialArray)rangedAgentManager.getAgent(atom.getType());
        IPotential[] potentials = potentialArray.getPotentials();
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
        potentialArray = getIntraPotentials((AtomTypeMolecule)atom.getType());
        potentials = potentialArray.getPotentials();
        for(int i=0; i<potentials.length; i++) {
            ((PotentialGroupNbr)potentials[i]).calculateRangeIndependent(atom,id,pc);
        }

        //cannot use AtomIterator field because of recursive call
        IAtomSet list = atom.getChildList();
        int size = list.getAtomCount();
        for (int i=0; i<size; i++) {
            calculate((IAtomLeaf)list.getAtom(i), id, pc);//recursive call
        }
	}
    
    /**
     * Performs given PotentialCalculation using potentials/neighbors associated
     * with the given atom (if any).  Then, if atom is not a leaf atom, iteration over
     * child atoms is performed and process is repeated (recursively) with each on down
     * the hierarchy until leaf atoms are reached.
     */
    //TODO make a "TerminalGroup" indicator in type that permits child atoms but indicates that no potentials apply directly to them
    protected void calculate(IAtomLeaf atom, IteratorDirective id, PotentialCalculation pc) {
        PotentialArray potentialArray = (PotentialArray)rangedAgentManager.getAgent(atom.getType());
        IPotential[] potentials = potentialArray.getPotentials();
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
    }
    
    private static final long serialVersionUID = 1L;
    private final AtomIteratorSinglet singletAtomIterator;
	private final AtomsetIteratorSinglet singletPairIterator;
    private int cellRange;
    protected final AtomsetIteratorPDT neighborIterator;
    private NeighborCriterion[] criteriaArray = new NeighborCriterion[0];
    
    public static class BoxAgentSiteManager implements BoxAgentSource {
        public BoxAgentSiteManager(int nCells, Space _space) {
            this.nCells = nCells;
            this.space = _space;
        }
        
        public Class getAgentClass() {
            return NeighborSiteManager.class;
        }
        
        public Object makeAgent(IBox box) {
            return new NeighborSiteManager(box,nCells, space);
        }
        
        public void releaseAgent(Object agent) {
        }
        
        private final int nCells;
        private final Space space;
    }
}
