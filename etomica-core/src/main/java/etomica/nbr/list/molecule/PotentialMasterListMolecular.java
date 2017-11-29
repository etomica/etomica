/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.list.molecule;

import etomica.box.Box;
import etomica.box.BoxAgentManager;
import etomica.box.BoxCellManager;
import etomica.molecule.*;
import etomica.molecule.iterator.MoleculeIteratorSinglet;
import etomica.nbr.cell.molecule.NeighborCellManagerMolecular;
import etomica.nbr.molecule.*;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.species.ISpecies;
import etomica.util.Arrays;
import etomica.util.Debug;

/**
 * PotentialMaster used to implement neighbor listing.  Instance of this
 * class is given as an argument to the Simulation constructor.
 * 
 * @author taitan
 *
 */
public class PotentialMasterListMolecular extends PotentialMasterNbrMolecular {

    /**
     * Default constructor uses range of 1.0.
     */
    public PotentialMasterListMolecular(Simulation sim, Space _space) {
        this(sim, 1.0, _space);
    }
    
    /**
     * Constructor specifying space and range for neighbor listing; uses null AtomPositionDefinition.
     */
    public PotentialMasterListMolecular(Simulation sim, double range, Space _space) {
        this(sim, range, (IMoleculePositionDefinition)null, _space);
    }
    
    /**
     * Constructs class using given position definition for all atom cell
     * assignments.
     * 
     * @param positionDefinition
     *            if null, specifies use of atom type's position definition
     */
    public PotentialMasterListMolecular(Simulation sim, double range, IMoleculePositionDefinition positionDefinition, Space _space) {
        this(sim, range, new BoxAgentSourceCellManagerListMolecular(sim, positionDefinition, _space), _space);
    }

    public PotentialMasterListMolecular(Simulation sim, double range, BoxAgentSourceCellManagerListMolecular boxAgentSource, Space _space) {
        this(sim, range, boxAgentSource, new BoxAgentManager<NeighborCellManagerMolecular>(boxAgentSource, NeighborCellManagerMolecular.class, sim), _space);
    }

    public PotentialMasterListMolecular(Simulation sim, double range, BoxAgentSourceCellManagerListMolecular boxAgentSource, BoxAgentManager<? extends BoxCellManager> agentManager, Space _space){
        this(sim, range, boxAgentSource, agentManager, new NeighborListAgentSourceMolecular(range, _space), _space);
    }

    public PotentialMasterListMolecular(Simulation sim, double range,
                                        BoxAgentSourceCellManagerListMolecular boxAgentSource,
                                        BoxAgentManager<? extends BoxCellManager> agentManager,
                                        NeighborListAgentSourceMolecular neighborListAgentSource, Space _space) {
        super(sim, boxAgentSource, agentManager);
        space = _space;
        this.neighborListAgentSource = neighborListAgentSource;
        neighborListAgentSource.setPotentialMaster(this);
        neighborListAgentManager = new BoxAgentManager<NeighborListManagerMolecular>(neighborListAgentSource, NeighborListManagerMolecular.class, sim);
        singletIterator = new MoleculeIteratorSinglet();
        moleculeSetSinglet = new MoleculeSetSinglet();
        moleculePair = new MoleculePair();
        cellRange = 2;
        allCriteria = new NeighborCriterionMolecular[0];

        neighborListAgentManager.setSimulation(sim);

        BoxAgentManager.AgentIterator<? extends BoxCellManager> iterator = boxAgentManager.makeIterator();
        iterator.reset();
        while (iterator.hasNext()) {
            NeighborCellManagerListMolecular cellManager = (NeighborCellManagerListMolecular)iterator.next();
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
        ((BoxAgentSourceCellManagerListMolecular)boxAgentSource).setRange(range);
        recomputeCriteriaRanges();
        
        BoxAgentManager.AgentIterator<? extends BoxCellManager> iterator = boxAgentManager.makeIterator();
        iterator.reset();
        while (iterator.hasNext()) {
            NeighborCellManagerMolecular cellManager = (NeighborCellManagerMolecular)iterator.next();
            cellManager.setPotentialRange(range);
        }

        neighborListAgentSource.setRange(newRange);
        
        BoxAgentManager.AgentIterator<NeighborListManagerMolecular> iteratorList = neighborListAgentManager.makeIterator();
        iteratorList.reset();
        while (iteratorList.hasNext()) {
            NeighborListManagerMolecular neighborListManager = iteratorList.next();
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
    protected void addRangedPotentialForSpecies(IPotentialMolecular potential, ISpecies[] species) {
        // we'll fix the neighbor range later in recomputeCriteriaRanges
        // 0 guarantees the simulation to be hosed if our range is less than the potential range
        // (since recomputeCriteriaRange will bail in that case)
        NeighborCriterionMolecular criterion;
        if (species.length == 2) {
            CriterionSimpleMolecular rangedCriterion = new CriterionSimpleMolecular(getSimulation(), space, potential.getRange(), 0.0);
            criterion = new CriterionSpeciesPair(rangedCriterion, species[0], species[1]);
         
            if (species[0] == species[1]) {
                criterion = new CriterionSimpleMolecular(getSimulation(), space, potential.getRange(), 0.0);
            }
        } else if(species.length == 1){
        	criterion = new CriterionSpecies(new CriterionAllMolecular(), species[0]);
        } else {
        	throw new RuntimeException("<PotentialMasterListMolecular> setRangedPotential!!!! ");
        }
        
        // add the criterion to all existing NeighborListManagers
        allCriteria = (NeighborCriterionMolecular[]) Arrays.addObject(allCriteria, criterion);
        
        for (int i=0; i<species.length; i++) {
            ((PotentialArrayMolecular)rangedAgentManager.getAgent(species[i])).setCriterion(potential, criterion);
        }
        if (potential.getRange() > maxPotentialRange) {
            maxPotentialRange = potential.getRange();
        }
        recomputeCriteriaRanges();

        BoxAgentManager.AgentIterator<NeighborListManagerMolecular> iterator = neighborListAgentManager.makeIterator();
        iterator.reset();
        while (iterator.hasNext()) {
            NeighborListManagerMolecular neighborListManager = iterator.next();
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
            PotentialArrayMolecular potentialArray = (PotentialArrayMolecular)rangedPotentialIterator.next();
            IPotential[] potentials = potentialArray.getPotentials();
            for (int i=0; i<potentials.length; i++) {
                if (potentials[i].getRange() > maxPotentialRange) {
                    maxPotentialRange = potentials[i].getRange();
                }
            }
        }
        recomputeCriteriaRanges();
        
        BoxAgentManager.AgentIterator<NeighborListManagerMolecular> iterator = neighborListAgentManager.makeIterator();
        iterator.reset();
        while (iterator.hasNext()) {
            NeighborListManagerMolecular neighborListManager = iterator.next();
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
        	PotentialArrayMolecular potentialArray = (PotentialArrayMolecular)rangedPotentialIterator.next();
            IPotential[] potentials = potentialArray.getPotentials();
            NeighborCriterionMolecular[] criteria = potentialArray.getCriteria();
            // this will double (or more) count criteria that apply to multiple atom types, but it won't hurt us
            for (int j=0; j<criteria.length; j++) {
                CriterionSimpleMolecular rangedCriterion = getRangedCriterion(criteria[j]);
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
    public NeighborCriterionMolecular getCriterion(IPotentialMolecular potential) {
        rangedPotentialIterator.reset();
        while (rangedPotentialIterator.hasNext()) {
            PotentialArrayMolecular potentialArray = (PotentialArrayMolecular)rangedPotentialIterator.next();
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
    private static CriterionSimpleMolecular getRangedCriterion(NeighborCriterionMolecular criterion) {
        if (criterion instanceof CriterionSimpleMolecular) {
            return (CriterionSimpleMolecular)criterion;
        }
        if (criterion instanceof CriterionAdapterMolecular) {
            return getRangedCriterion(((CriterionAdapterMolecular)criterion).getWrappedCriterion());
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
    public void setCriterion(IPotentialMolecular potential, NeighborCriterionMolecular criterion) {
        NeighborCriterionMolecular oldCriterion = getCriterion(potential);
        if (oldCriterion != null) {
            // remove the criterion to all existing NeighborListManagers
            allCriteria = (NeighborCriterionMolecular[]) Arrays.removeObject(allCriteria, criterion);
        }
        rangedPotentialIterator.reset();
        boolean success = false;
        while (rangedPotentialIterator.hasNext()) {
            PotentialArrayMolecular potentialArray = (PotentialArrayMolecular)rangedPotentialIterator.next();
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
            allCriteria = (NeighborCriterionMolecular[]) Arrays.addObject(allCriteria, criterion);
        	return;
        }
        throw new IllegalArgumentException("Potential "+potential+" is not associated with this PotentialMasterList");
    }
    
    public void removePotential(IPotentialMolecular potential) {
        rangedPotentialIterator.reset();
        while (rangedPotentialIterator.hasNext()) {
            PotentialArrayMolecular potentialArray = (PotentialArrayMolecular)rangedPotentialIterator.next();
            IPotential[] potentials = potentialArray.getPotentials();
            for (int j=0; j<potentials.length; j++) {
                if (potentials[j] == potential) {
                    // found it!
                    // remove the criterion from our list
                    allCriteria = (NeighborCriterionMolecular[]) Arrays.removeObject(allCriteria, potentialArray.getCriteria()[j]);
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

        BoxAgentManager.AgentIterator<NeighborListManagerMolecular> iterator = neighborListAgentManager.makeIterator();
        iterator.reset();
        while (iterator.hasNext()) {
            NeighborListManagerMolecular neighborListManager = iterator.next();
            neighborListManager.updateLists();
        }
    }
    
    public NeighborCriterionMolecular[] getNeighborCriteria() {
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
        IMolecule targetMolecule = id.getTargetMolecule();
        NeighborListManagerMolecular neighborManager = neighborListAgentManager.getAgent(box);
        if (targetMolecule == null) {
            if (Debug.ON && id.direction() != IteratorDirective.Direction.UP) {
                throw new IllegalArgumentException("When there is no target, iterator directive must be up");
            }
            // invoke setBox on all potentials
            for (int i=0; i<allPotentials.length; i++) {
                allPotentials[i].setBox(box);
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
      
        	PotentialArrayMolecular potentialArray = getIntraPotentials(targetMolecule.getType());
        	IPotential[] potentials = potentialArray.getPotentials();
        	for(int i=0; i<potentials.length; i++) {
        		potentials[i].setBox(box);
        	}

        	calculate(targetMolecule, id.direction(), pc, neighborManager);
            
        }
        if(lrcMaster != null) {
            lrcMaster.calculate(box, id, pc);
        }
    }
    
    protected void calculate(IMolecule molecule, IteratorDirective.Direction direction, PotentialCalculation pc, NeighborListManagerMolecular neighborManager) {
    	//System.out.println("here");
    	singletIterator.setMolecule(molecule);
        PotentialArrayMolecular potentialArray = (PotentialArrayMolecular)rangedAgentManager.getAgent(molecule.getType());
        IPotential[] potentials = potentialArray.getPotentials();
        for(int i=0; i<potentials.length; i++) {
            switch (potentials[i].nBody()) {
            case 1:
                boolean[] potential1BodyArray = neighborManager.getPotential1BodyList(molecule).getInteractingList();
                if (potential1BodyArray[i]) {
                    moleculeSetSinglet.atom = molecule;
                    ((PotentialCalculationMolecular)pc).doCalculation(moleculeSetSinglet, (IPotentialMolecular)potentials[i]);
                }
                break;
            case 2:
                if (direction != IteratorDirective.Direction.DOWN) {
                    IMoleculeList list = neighborManager.getUpList(molecule)[i];
                    int nNeighbors = list.getMoleculeCount();
                    moleculePair.atom0 = molecule;
                    for (int j=0; j<nNeighbors; j++) {
                        moleculePair.atom1 = list.getMolecule(j);
                        ((PotentialCalculationMolecular)pc).doCalculation(moleculePair, (IPotentialMolecular)potentials[i]);
                    }
                }
                if (direction != IteratorDirective.Direction.UP) {
                    IMoleculeList list = neighborManager.getDownList(molecule)[i];
                    int nNeighbors = list.getMoleculeCount();
                    moleculePair.atom1 = molecule;
                    for (int j=0; j<nNeighbors; j++) {
                        moleculePair.atom0 = list.getMolecule(j);
                        ((PotentialCalculationMolecular)pc).doCalculation(moleculePair, (IPotentialMolecular)potentials[i]);
                    }
                }
                break;//switch
            case Integer.MAX_VALUE: //N-body
                // do the calculation considering the current Atom as the 
                // "central" Atom.
            	if(moleculeArrayList==null){
            		moleculeArrayList = new MoleculeArrayList();
            	}
                doNBodyStuff(molecule, pc, i, (IPotentialMolecular)potentials[i], neighborManager);
                if (direction != IteratorDirective.Direction.UP) {
                    // must have a target and be doing "both"
                    // we have to do the calculation considering each of the 
                    // target's neighbors
                    IMoleculeList list = neighborManager.getUpList(molecule)[i];
                    for (int j=0; j<list.getMoleculeCount(); j++) {
                        IMolecule otherMolecule = list.getMolecule(j);
                        doNBodyStuff(otherMolecule, pc, i, (IPotentialMolecular)potentials[i], neighborManager);
                    }
                    list = neighborManager.getDownList(molecule)[i];
                    for (int j=0; j<list.getMoleculeCount(); j++) {
                        IMolecule otherMolecule = list.getMolecule(j);
                        doNBodyStuff(otherMolecule, pc, i, (IPotentialMolecular)potentials[i], neighborManager);
                    }
                }
                
            }//end of switch
        }//end of for
    }
    
    /**
     * Invokes the PotentialCalculationMolecular for the given Molecule with its up and down
     * neighbors as a single MoleculeSet.
     */
    protected void doNBodyStuff(IMolecule molecule, PotentialCalculation pc, int potentialIndex, 
            IPotentialMolecular potential, NeighborListManagerMolecular neighborManager) {
        moleculeArrayList.add(molecule);
        IMoleculeList[] list = neighborManager.getUpList(molecule);
        if (potentialIndex < list.length) {
            moleculeArrayList.addAll(list[potentialIndex]);
        }
        list = neighborManager.getDownList(molecule);
        if (potentialIndex < list.length) {
            moleculeArrayList.addAll(list[potentialIndex]);
        }
        ((PotentialCalculationMolecular)pc).doCalculation(moleculeArrayList, potential);
        moleculeArrayList.clear();
    }

    public NeighborListManagerMolecular getNeighborManager(Box box) {
        // we didn't have the simulation when we made the agent manager.
        // setting the simulation after the first time is a quick return
        return neighborListAgentManager.getAgent(box);
    }

    
    public NeighborCellManagerMolecular getNbrCellManager(Box box) {
        return (NeighborCellManagerMolecular)boxAgentManager.getAgent(box);
    }

    public void setCellRange(int newCellRange) {
        cellRange = newCellRange;

        BoxAgentManager.AgentIterator<? extends BoxCellManager> iterator = boxAgentManager.makeIterator();
        iterator.reset();
        while (iterator.hasNext()) {
            NeighborCellManagerMolecular cellManager = (NeighborCellManagerMolecular)iterator.next();
            cellManager.setCellRange(cellRange);
        }
    }

    public int getCellRange() {
        return cellRange;
    }

    protected final Space space;
    private final MoleculeIteratorSinglet singletIterator;
    protected final MoleculeSetSinglet moleculeSetSinglet;
    protected final MoleculePair moleculePair;
    protected final NeighborListAgentSourceMolecular neighborListAgentSource;
    protected final BoxAgentManager<NeighborListManagerMolecular> neighborListAgentManager;
    private int cellRange;
    protected double range;
    private double maxPotentialRange = 0;
    private double safetyFactor = 0.4;
    protected NeighborCriterionMolecular[] allCriteria;
    
    // things needed for N-body potentials
    private MoleculeArrayList moleculeArrayList;
    
    public static class NeighborListAgentSourceMolecular implements BoxAgentManager.BoxAgentSource<NeighborListManagerMolecular>{
        public NeighborListAgentSourceMolecular(double range, Space space) {
            
            this.range = range;
            this.space = space;
        }
        
        public void setRange(double newRange) {
            range = newRange;
        }
        
        public void setPotentialMaster(PotentialMasterListMolecular p){
            potentialMaster = p;
        }

        public NeighborListManagerMolecular makeAgent(Box box) {
            return new NeighborListManagerMolecular(potentialMaster, range, box, space);
        }
        
        public void releaseAgent(NeighborListManagerMolecular object) {
            object.dispose();
        }
        
        protected PotentialMasterListMolecular potentialMaster;
        protected double range;
        protected final Space space;
    }


}
