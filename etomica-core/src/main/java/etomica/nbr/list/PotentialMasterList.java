/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.list;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.IAtomType;
import etomica.box.Box;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.IPotential;
import etomica.api.IPotentialAtomic;
import etomica.simulation.Simulation;
import etomica.api.ISpecies;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomPair;
import etomica.atom.AtomSetSinglet;
import etomica.atom.IAtomPositionDefinition;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.atom.iterator.IteratorDirective;
import etomica.box.BoxAgentManager;
import etomica.box.BoxCellManager;
import etomica.nbr.CriterionAdapter;
import etomica.nbr.CriterionAll;
import etomica.nbr.CriterionInterMolecular;
import etomica.nbr.CriterionSimple;
import etomica.nbr.CriterionType;
import etomica.nbr.CriterionTypePair;
import etomica.nbr.CriterionTypesMulti;
import etomica.nbr.NeighborCriterion;
import etomica.nbr.PotentialGroupNbr;
import etomica.nbr.PotentialMasterNbr;
import etomica.nbr.cell.NeighborCellManager;
import etomica.potential.PotentialArray;
import etomica.potential.PotentialCalculation;
import etomica.potential.PotentialGroup;
import etomica.space.Space;
import etomica.util.Arrays;
import etomica.util.Debug;

/**
 * PotentialMaster used to implement neighbor listing.  Instance of this
 * class is given as an argument to the Simulation constructor.
 */
public class PotentialMasterList extends PotentialMasterNbr {

    /**
     * Default constructor uses range of 1.0.
     */
    public PotentialMasterList(Simulation sim, Space _space) {
        this(sim,1.0, _space);
    }
    
    /**
     * Constructor specifying space and range for neighbor listing; uses null AtomPositionDefinition.
     */
    public PotentialMasterList(Simulation sim, double range, Space _space) {
        this(sim, range, (IAtomPositionDefinition)null, _space);
    }
    
    /**
     * Constructs class using given position definition for all atom cell
     * assignments.
     * 
     * @param positionDefinition
     *            if null, specifies use of atom type's position definition
     */
    public PotentialMasterList(Simulation sim, double range, IAtomPositionDefinition positionDefinition, Space _space) {
        this(sim, range, new BoxAgentSourceCellManagerList(sim, positionDefinition, _space), _space);
    }

    public PotentialMasterList(Simulation sim, double range, BoxAgentSourceCellManagerList boxAgentSource, Space _space) {
        this(sim, range, boxAgentSource, new BoxAgentManager<NeighborCellManager>(boxAgentSource, NeighborCellManager.class), _space);
    }

    public PotentialMasterList(Simulation sim, double range, BoxAgentSourceCellManagerList boxAgentSource, BoxAgentManager<? extends BoxCellManager> agentManager, Space _space){
        this(sim, range, boxAgentSource, agentManager, new NeighborListAgentSource(range, _space), _space);
    }

    public PotentialMasterList(Simulation sim, double range,
                               BoxAgentSourceCellManagerList boxAgentSource,
                               BoxAgentManager<? extends BoxCellManager> agentManager,
                               NeighborListAgentSource neighborListAgentSource, Space _space) {
        super(sim, boxAgentSource, agentManager);
        space = _space;
        this.neighborListAgentSource = neighborListAgentSource;
        neighborListAgentSource.setPotentialMaster(this);
        neighborListAgentManager = new BoxAgentManager<NeighborListManager>(neighborListAgentSource, NeighborListManager.class);
        singletIterator = new AtomIteratorSinglet();
        atomSetSinglet = new AtomSetSinglet();
        atomPair = new AtomPair();
        cellRange = 2;
        allCriteria = new NeighborCriterion[0];

        neighborListAgentManager.setSimulation(sim);

        BoxAgentManager.AgentIterator<? extends BoxCellManager> iterator = boxAgentManager.makeIterator();
        iterator.reset();
        while (iterator.hasNext()) {
            NeighborCellManagerList cellManager = (NeighborCellManagerList)iterator.next();
            cellManager.setPotentialMaster(this);
        }

        boxAgentSource.setPotentialMaster(this);

        // setRange last.  that should always be OK since anyone can call
        // setRange later. if we call it early, member fields won't exist,
        // if we just set range, our stuff associated with existing boxes
        // won't get initialized
        setRange(range);
    }
    
    /**
     * Sets the range that determines how far to look for neighbors.  This 
     * range must be greater than the range of the longest-range potential.
     */
    public void setRange(double newRange) {
        if (newRange <= 0) {
            throw new IllegalArgumentException("Range must be greater than 0");
        }
        range = newRange;
        ((BoxAgentSourceCellManagerList)boxAgentSource).setRange(range);
        recomputeCriteriaRanges();
        
        BoxAgentManager.AgentIterator<? extends BoxCellManager> iterator = boxAgentManager.makeIterator();
        iterator.reset();
        while (iterator.hasNext()) {
            NeighborCellManager cellManager = (NeighborCellManager)iterator.next();
            cellManager.setPotentialRange(range);
        }

        neighborListAgentSource.setRange(newRange);
        
        BoxAgentManager.AgentIterator<NeighborListManager> iteratorList = neighborListAgentManager.makeIterator();
        iteratorList.reset();
        while (iteratorList.hasNext()) {
            NeighborListManager neighborListManager = iteratorList.next();
            neighborListManager.setRange(range);
        }
    }
    
    /**
     * Returns the range that determines how far to look for neighbors.
     */
    public double getRange() {
        return range;
    }
    
    /**
     * Returns the maximum range of any potential held by this potential master
     */
    public double getMaxPotentialRange() {
        return maxPotentialRange;
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
        if (newSafetyFactor <= 0 || newSafetyFactor > 0.5) {
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
     * Add the given potential to be used for the given atom types and the
     * given criterion.  If multiple types are given, then the potential will
     * be used for any set of atoms containing only the atom types listed;
     * if types A and B are passed, then the potential will be used for
     * A-A, A-B, B-A and B-B.  The criterion must handle any additional
     * filtering.
     * 
     * The given potential will not be held by a PotentialGroup.
     */
    public void addPotentialNbrList(IPotentialAtomic potential, IAtomType[] atomTypes, NeighborCriterion criterion) {
        if (potential.getRange() == Double.POSITIVE_INFINITY) {
            throw new RuntimeException("not the method you wanted to call");
        }
        for (int i=0; i<atomTypes.length; i++) {
            addRangedPotential(potential, atomTypes[i]);
        }

        // add the criterion to all existing NeighborListManagers
        allCriteria = (NeighborCriterion[]) Arrays.addObject(allCriteria, criterion);
        
        for (int i=0; i<atomTypes.length; i++) {
            ((PotentialArray)rangedAgentManager.getAgent(atomTypes[i])).setCriterion(potential, criterion);
        }
        if (potential.getRange() > maxPotentialRange) {
            maxPotentialRange = potential.getRange();
        }
        recomputeCriteriaRanges();

        BoxAgentManager.AgentIterator<NeighborListManager> iterator = neighborListAgentManager.makeIterator();
        iterator.reset();
        while (iterator.hasNext()) {
            NeighborListManager neighborListManager = iterator.next();
            neighborListManager.updateLists();
        }
    }

    /**
     * Adds the potential as a ranged potential that applies to the given 
     * AtomTypes.  This method creates a criterion for the potential and 
     * notifies the NeighborListManager of its existence.
     */
    protected void addRangedPotentialForTypes(IPotentialAtomic potential, IAtomType[] atomType) {
        // we'll fix the neighbor range later in recomputeCriteriaRanges
        // 0 guarantees the simulation to be hosed if our range is less than the potential range
        // (since recomputeCriteriaRange will bail in that case)
        NeighborCriterion criterion;
        if (atomType.length == 2) {
            NeighborCriterion rangedCriterion;
            if (potential.getRange() < Double.POSITIVE_INFINITY) {
                rangedCriterion = new CriterionSimple(getSimulation(), space, potential.getRange(), 0.0);
            }
            else {
                rangedCriterion = new CriterionAll();
            }
            criterion = new CriterionTypePair(rangedCriterion, atomType[0], atomType[1]);
            ISpecies moleculeType0 = atomType[0].getSpecies();
            ISpecies moleculeType1 = atomType[1].getSpecies();
            if (moleculeType0 == moleculeType1) {
                criterion = new CriterionInterMolecular(criterion);
            }
        }
        else if (atomType.length == 1) {
            criterion = new CriterionType(new CriterionAll(), atomType[0]);
        }
        else {
            criterion = new CriterionTypesMulti(new CriterionAll(), atomType);
        }

        // add the criterion to all existing NeighborListManagers
        allCriteria = (NeighborCriterion[]) Arrays.addObject(allCriteria, criterion);
        
        for (int i=0; i<atomType.length; i++) {
            ((PotentialArray)rangedAgentManager.getAgent(atomType[i])).setCriterion(potential, criterion);
        }
        if (potential.getRange() > maxPotentialRange && potential.getRange() < Double.POSITIVE_INFINITY) {
            maxPotentialRange = potential.getRange();
        }
        recomputeCriteriaRanges();

        BoxAgentManager.AgentIterator<NeighborListManager> iterator = neighborListAgentManager.makeIterator();
        iterator.reset();
        while (iterator.hasNext()) {
            NeighborListManager neighborListManager = iterator.next();
            neighborListManager.updateLists();
        }
    }
    
    /**
     * Recomputes the maximum potential range (which might change without this
     * class receiving notification) and readjust cell lists
     */
    public void reset() {
        rangedPotentialIterator.reset();
        maxPotentialRange = 0;
        while (rangedPotentialIterator.hasNext()) {
            PotentialArray potentialArray = (PotentialArray)rangedPotentialIterator.next();
            IPotential[] potentials = potentialArray.getPotentials();
            for (int i=0; i<potentials.length; i++) {
                if (potentials[i].getRange() > maxPotentialRange) {
                    maxPotentialRange = potentials[i].getRange();
                }
            }
        }
        recomputeCriteriaRanges();
        
        BoxAgentManager.AgentIterator<NeighborListManager> iterator = neighborListAgentManager.makeIterator();
        iterator.reset();
        while (iterator.hasNext()) {
            NeighborListManager neighborListManager = iterator.next();
            neighborListManager.reset();
        }
    }
    
    /**
     * Recomputes the range for all criterion based on our own range.  The 
     * range for each criterion is set so that they all have an equal max
     * displacement.  Our nominal neighbor range is used for the criterion
     * with the longest potential range.
     */
    public void recomputeCriteriaRanges() {
        double maxDisplacement = (getRange() - maxPotentialRange) * safetyFactor;
        if (maxDisplacement < 0) {
            // someone probably added a long ranged potential and hasn't updated the PotentialMaster's range
            // when they do, we'll get called again.
            // if they don't, the simulation will probably crash
            return;
        }
        rangedPotentialIterator.reset();
        while (rangedPotentialIterator.hasNext()) {
            PotentialArray potentialArray = (PotentialArray)rangedPotentialIterator.next();
            IPotential[] potentials = potentialArray.getPotentials();
            NeighborCriterion[] criteria = potentialArray.getCriteria();
            // this will double (or more) count criteria that apply to multiple atom types, but it won't hurt us
            for (int j=0; j<criteria.length; j++) {
                CriterionSimple rangedCriterion = getRangedCriterion(criteria[j]);
                if (rangedCriterion != null) {
                    double newRange = maxDisplacement/safetyFactor + potentials[j].getRange();
                    rangedCriterion.setNeighborRange(newRange);
                    rangedCriterion.setInteractionRange(potentials[j].getRange());
                    rangedCriterion.setSafetyFactor(safetyFactor);
                }                    
            }
        }
    }
    
    /**
     * Returns the criterion used by to determine what atoms interact with the
     * given potential.
     */
    public NeighborCriterion getCriterion(IPotentialAtomic potential) {
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
    public void setCriterion(IPotentialAtomic potential, NeighborCriterion criterion) {
        NeighborCriterion oldCriterion = getCriterion(potential);
        if (oldCriterion != null) {
            // remove the criterion to all existing NeighborListManagers
            allCriteria = (NeighborCriterion[]) Arrays.removeObject(allCriteria, oldCriterion);
        }
        rangedPotentialIterator.reset();
        boolean success = false;
        while (rangedPotentialIterator.hasNext()) {
            PotentialArray potentialArray = (PotentialArray)rangedPotentialIterator.next();
            IPotential[] potentials = potentialArray.getPotentials();
            for (int j=0; j<potentials.length; j++) {
                if (potentials[j] == potential) {
                    success = true;
                    potentialArray.setCriterion(potential, criterion);
                    break;
                }
            }
        }
        if (success) {
            // add the criterion to all existing NeighborListManagers
            allCriteria = (NeighborCriterion[]) Arrays.addObject(allCriteria, criterion);
        	return;
        }
        throw new IllegalArgumentException("Potential "+potential+" is not associated with this PotentialMasterList");
    }
    
    public void removePotential(IPotentialAtomic potential) {
        rangedPotentialIterator.reset();
        while (rangedPotentialIterator.hasNext()) {
            PotentialArray potentialArray = (PotentialArray)rangedPotentialIterator.next();
            IPotential[] potentials = potentialArray.getPotentials();
            for (int j=0; j<potentials.length; j++) {
                if (potentials[j] == potential) {
                    // found it!
                    // remove the criterion from our list
                    allCriteria = (NeighborCriterion[]) Arrays.removeObject(allCriteria, potentialArray.getCriteria()[j]);
                    break;
                }
            }
        }

        super.removePotential(potential);
        
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

        BoxAgentManager.AgentIterator<NeighborListManager> iterator = neighborListAgentManager.makeIterator();
        iterator.reset();
        while (iterator.hasNext()) {
            NeighborListManager neighborListManager = iterator.next();
            neighborListManager.updateLists();
        }
    }
    
    public NeighborCriterion[] getNeighborCriteria() {
        return allCriteria;
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
    public void calculate(Box box, IteratorDirective id, PotentialCalculation pc) {
        if(!enabled) return;
        IAtom targetAtom = id.getTargetAtom();
        IMolecule targetMolecule = id.getTargetMolecule();
        NeighborListManager neighborManager = neighborListAgentManager.getAgent(box);
        if (targetAtom == null && targetMolecule == null) {
            if (Debug.ON && id.direction() != IteratorDirective.Direction.UP) {
                throw new IllegalArgumentException("When there is no target, iterator directive must be up");
            }
            // invoke setBox on all potentials
            for (int i=0; i<allPotentials.length; i++) {
                allPotentials[i].setBox(box);
                if(allPotentials[i].nBody() == 0){
                	((PotentialGroup)allPotentials[i]).calculate(new MoleculeIterator0(), id.direction(), null, pc);
                }
            }

            //no target atoms specified
            //call calculate with each SpeciesAgent
            IMoleculeList list = box.getMoleculeList();
            int size = list.getMoleculeCount();
            for (int i=0; i<size; i++) {
                calculate(list.getMolecule(i), id.direction(), pc, neighborManager);//call calculate with the SpeciesAgent
            }
        }
        else {
            if (targetAtom != null) {
                PotentialArray potentialArray = (PotentialArray)rangedAgentManager.getAgent(targetAtom.getType());
                IPotential[] potentials = potentialArray.getPotentials();
                for(int i=0; i<potentials.length; i++) {
                    potentials[i].setBox(box);
                }
                
                //first walk up the tree looking for 1-body range-independent potentials that apply to parents
                IMolecule parentAtom = targetAtom.getParentGroup();
                potentialArray = getIntraPotentials(parentAtom.getType());
                potentials = potentialArray.getPotentials();
                for(int i=0; i<potentials.length; i++) {
                    potentials[i].setBox(box);
                    ((PotentialGroupNbr)potentials[i]).calculateRangeIndependent(parentAtom,id.direction(), targetAtom,pc);
                }
                calculate(targetAtom, id.direction(), pc, neighborManager);
            }
            else {
                PotentialArray potentialArray = getIntraPotentials(targetMolecule.getType());
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
	
    /**
     * Performs given PotentialCalculation using potentials/neighbors associated
     * with the given atom (if any).  Then, if atom is not a leaf atom, iteration over
     * child atoms is performed and process is repeated (recursively) with each on down
     * the hierarchy until leaf atoms are reached.
     */
    protected void calculate(IMolecule molecule, IteratorDirective.Direction direction, PotentialCalculation pc, NeighborListManager neighborManager) {
        PotentialArray potentialArray = getIntraPotentials(molecule.getType());
        IPotential[] potentials = potentialArray.getPotentials();
        for(int i=0; i<potentials.length; i++) {
            ((PotentialGroupNbr)potentials[i]).calculateRangeIndependent(molecule,direction, null, pc);
        }
        
        //cannot use AtomIterator field because of recursive call
        IAtomList list = molecule.getChildList();
        int size = list.getAtomCount();
        for (int i=0; i<size; i++) {
            calculate(list.getAtom(i), direction, pc, neighborManager);//recursive call
        }
    }
    
    protected void calculate(IAtom atom, IteratorDirective.Direction direction, PotentialCalculation pc, NeighborListManager neighborManager) {
        singletIterator.setAtom(atom);
        PotentialArray potentialArray = (PotentialArray)rangedAgentManager.getAgent(atom.getType());
        IPotential[] potentials = potentialArray.getPotentials();
        for(int i=0; i<potentials.length; i++) {
            switch (potentials[i].nBody()) {
            case 1:
                boolean[] potential1BodyArray = neighborManager.getPotential1BodyList(atom).getInteractingList();
                if (potential1BodyArray[i]) {
                    atomSetSinglet.atom = atom;
                    pc.doCalculation(atomSetSinglet, (IPotentialAtomic)potentials[i]);
                }
                break;
            case 2:
                if (direction != IteratorDirective.Direction.DOWN) {
                    IAtomList list = neighborManager.getUpList(atom)[i];
                    int nNeighbors = list.getAtomCount();
                    atomPair.atom0 = atom;
                    for (int j=0; j<nNeighbors; j++) {
                        atomPair.atom1 = list.getAtom(j);
                        pc.doCalculation(atomPair, (IPotentialAtomic)potentials[i]);
                    }
                }
                if (direction != IteratorDirective.Direction.UP) {
                    IAtomList list = neighborManager.getDownList(atom)[i];
                    int nNeighbors = list.getAtomCount();
                    atomPair.atom1 = atom;
                    for (int j=0; j<nNeighbors; j++) {
                        atomPair.atom0 = list.getAtom(j);
                        pc.doCalculation(atomPair, (IPotentialAtomic)potentials[i]);
                    }
                }
                break;//switch
            case Integer.MAX_VALUE: //N-body
                // do the calculation considering the current Atom as the 
                // "central" Atom.
            	if(atomArrayList==null){
            		atomArrayList = new AtomArrayList();
            	}
                doNBodyStuff(atom, pc, i, (IPotentialAtomic)potentials[i], neighborManager);
                if (direction != IteratorDirective.Direction.UP) {
                    // must have a target and be doing "both"
                    // we have to do the calculation considering each of the 
                    // target's neighbors
                    IAtomList list = neighborManager.getUpList(atom)[i];
                    for (int j=0; j<list.getAtomCount(); j++) {
                        IAtom otherAtom = list.getAtom(j);
                        doNBodyStuff(otherAtom, pc, i, (IPotentialAtomic)potentials[i], neighborManager);
                    }
                    list = neighborManager.getDownList(atom)[i];
                    for (int j=0; j<list.getAtomCount(); j++) {
                        IAtom otherAtom = list.getAtom(j);
                        doNBodyStuff(otherAtom, pc, i, (IPotentialAtomic)potentials[i], neighborManager);
                    }
                }
                
            }//end of switch
        }//end of for
    }
    
    /**
     * Invokes the PotentialCalculation for the given Atom with its up and down
     * neighbors as a single AtomSet.
     */
    protected void doNBodyStuff(IAtom atom, PotentialCalculation pc, int potentialIndex, 
            IPotentialAtomic potential, NeighborListManager neighborManager) {
        atomArrayList.add(atom);
        IAtomList[] list = neighborManager.getUpList(atom);
        if (potentialIndex < list.length) {
            atomArrayList.addAll(list[potentialIndex]);
        }
        list = neighborManager.getDownList(atom);
        if (potentialIndex < list.length) {
            atomArrayList.addAll(list[potentialIndex]);
        }
        pc.doCalculation(atomArrayList, potential);
        atomArrayList.clear();
    }

    public NeighborListManager getNeighborManager(Box box) {
        // we didn't have the simulation when we made the agent manager.
        // setting the simulation after the first time is a quick return
        return neighborListAgentManager.getAgent(box);
    }

    
    public NeighborCellManager getNbrCellManager(Box box) {
        return (NeighborCellManager)boxAgentManager.getAgent(box);
    }

    public void setCellRange(int newCellRange) {
        cellRange = newCellRange;

        BoxAgentManager.AgentIterator<? extends BoxCellManager> iterator = boxAgentManager.makeIterator();
        iterator.reset();
        while (iterator.hasNext()) {
            NeighborCellManager cellManager = (NeighborCellManager)iterator.next();
            cellManager.setCellRange(cellRange);
        }
    }

    public int getCellRange() {
        return cellRange;
    }

    protected final Space space;
    private final AtomIteratorSinglet singletIterator;
    protected final AtomSetSinglet atomSetSinglet;
    protected final AtomPair atomPair;
    protected final NeighborListAgentSource neighborListAgentSource;
    protected final BoxAgentManager<NeighborListManager> neighborListAgentManager;
    private int cellRange;
    protected double range;
    private double maxPotentialRange = 0;
    private double safetyFactor = 0.4;
    protected NeighborCriterion[] allCriteria;
    
    // things needed for N-body potentials
    private AtomArrayList atomArrayList;
    
    public static class NeighborListAgentSource implements BoxAgentManager.BoxAgentSource<NeighborListManager> {
        public NeighborListAgentSource(double range, Space space) {
            
            this.range = range;
            this.space = space;
        }
        
        public void setRange(double newRange) {
            range = newRange;
        }
        
        public void setPotentialMaster(PotentialMasterList p){
            potentialMaster = p;
        }
        
        public NeighborListManager makeAgent(Box box) {
            return new NeighborListManager(potentialMaster, range, box, space);
        }
        
        public void releaseAgent(NeighborListManager object) {
            object.dispose();
        }
        
        protected PotentialMasterList potentialMaster;
        protected double range;
        protected final Space space;
    }
}
