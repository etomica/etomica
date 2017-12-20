/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.site;

import etomica.atom.AtomSetSinglet;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.iterator.AtomsetIteratorPDT;
import etomica.box.Box;
import etomica.box.BoxAgentManager;
import etomica.box.BoxAgentManager.BoxAgentSource;
import etomica.box.BoxCellManager;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.nbr.*;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.species.ISpecies;
import etomica.util.Arrays;
import etomica.util.Debug;

public class PotentialMasterSite extends PotentialMasterNbr {

    protected final AtomSetSinglet atomSetSinglet;
    protected final AtomsetIteratorPDT neighborIterator;
    private int cellRange;
    private NeighborCriterion[] criteriaArray = new NeighborCriterion[0];
    
	/**
	 * Invokes superclass constructor, specifying IteratorFactoryCell
     * for generating molecule iterators.  Sets default nCells of 10 and
     * position definition to null, so that atom type's definition is used
     * to assign cells.
     */
	public PotentialMasterSite(Simulation sim, int nCells, Space _space) {
        this(sim, new BoxAgentSiteSource(nCells, _space), _space);
    }

    public PotentialMasterSite(Simulation sim,
    		                   BoxAgentSource<BoxCellManager> boxAgentSource, Space _space) {
        this(sim, boxAgentSource, new BoxAgentManager<BoxCellManager>(boxAgentSource, BoxCellManager.class, sim), _space);
    }
    
    public PotentialMasterSite(Simulation sim, BoxAgentSource<BoxCellManager> boxAgentSource,
                               BoxAgentManager<BoxCellManager> agentManager, Space _space) {
        this(sim, boxAgentSource, agentManager, new Api1ASite(_space.D(),agentManager));
    }
    
    protected PotentialMasterSite(Simulation sim, BoxAgentSource<? extends BoxCellManager> boxAgentSource,
                                  BoxAgentManager<? extends BoxCellManager> agentManager, AtomsetIteratorPDT neighborIterator) {
        super(sim, boxAgentSource, agentManager);
        atomSetSinglet = new AtomSetSinglet();
        this.neighborIterator = neighborIterator;
	}
    
    /**
     * @return Returns the cellRange.
     */
    public int getCellRange() {
        return cellRange;
    }

    /**
     * @param newCellRange The cellRange to set.
     */
    public void setCellRange(int newCellRange) {
        cellRange = newCellRange;
    }

    protected void addRangedPotentialForTypes(IPotentialAtomic potential, AtomType[] atomType) {
        NeighborCriterion criterion;
        if (atomType.length == 2) {
            criterion = new CriterionTypePair(new CriterionAll(), atomType[0], atomType[1]);
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
        for (int i=0; i<atomType.length; i++) {
            ((PotentialArray)rangedAgentManager.getAgent(atomType[i])).setCriterion(potential, criterion);
        }
        criteriaArray = (NeighborCriterion[]) Arrays.addObject(criteriaArray, criterion);
    }

    /**
     * Returns the criterion used by to determine what atoms interact with the
     * given potential.
     */
    public NeighborCriterion getCriterion(IPotentialAtomic potential) {
        for (PotentialArray potentialArray : this.rangedAgentManager.getAgents().values()) {
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
    public void setCriterion(IPotentialAtomic potential, NeighborCriterion criterion) {
        for (PotentialArray potentialArray : this.rangedAgentManager.getAgents().values()) {
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
    public void calculate(Box box, IteratorDirective id, PotentialCalculation pc) {
        if (!enabled)
            return;
        for (int i=0; i<criteriaArray.length; i++) {
            criteriaArray[i].setBox(box);
        }
        IAtom targetAtom = id.getTargetAtom();
        IMolecule targetMolecule = id.getTargetMolecule();
        neighborIterator.setBox(box);
        if (targetAtom == null && targetMolecule == null) {
            if (Debug.ON && id.direction() != IteratorDirective.Direction.UP) {
                throw new IllegalArgumentException("When there is no target, iterator directive must be up");
            }
            neighborIterator.setDirection(IteratorDirective.Direction.UP);
            // invoke setBox on all potentials
            for (int i=0; i<allPotentials.length; i++) {
                allPotentials[i].setBox(box);
            }

            //no target atoms specified
            //call calculate with each molecule
            for (int j=0; j<simulation.getSpeciesCount(); j++) {
                IMoleculeList moleculeList = box.getMoleculeList(simulation.getSpecies(j));
                int size = moleculeList.getMoleculeCount();
                PotentialArray intraPotentialArray = getIntraPotentials(simulation.getSpecies(j));
                final IPotential[] intraPotentials = intraPotentialArray.getPotentials();
                for (int i=0; i<size; i++) {
                    IMolecule molecule = moleculeList.getMolecule(i);

                    IAtomList atomList = molecule.getChildList();
                    int numChildren = atomList.getAtomCount();
                    for (int k=0; k<numChildren; k++) {
                        calculate(atomList.getAtom(k), pc);
                    }

                    for(int k=0; k<intraPotentials.length; k++) {
                        ((PotentialGroupNbr)intraPotentials[k]).calculateRangeIndependent(molecule,id.direction(), null, pc);
                    }
                }
            }
        }
        else {
            // one target atom
            neighborIterator.setDirection(id.direction());
            if (targetAtom != null) {
                PotentialArray potentialArray = (PotentialArray)rangedAgentManager.getAgent(targetAtom.getType());
                IPotential[] potentials = potentialArray.getPotentials();
                for(int i=0; i<potentials.length; i++) {
                    potentials[i].setBox(box);
                }

                //walk up the tree looking for 1-body range-independent potentials that apply to parents
                IMolecule parentMolecule = targetAtom.getParentGroup();
                potentialArray = getIntraPotentials(parentMolecule.getType());
                potentials = potentialArray.getPotentials();
                for(int i=0; i<potentials.length; i++) {
                    potentials[i].setBox(box);
                    ((PotentialGroupNbr)potentials[i]).calculateRangeIndependent(parentMolecule,id.direction(), targetAtom, pc);
                }
                calculate(targetAtom, pc);
            }
            else {
                for (int i=0; i<allPotentials.length; i++) {
                    allPotentials[i].setBox(box);
                }

                IAtomList atomList = (targetMolecule).getChildList();
                int numChildren = atomList.getAtomCount();
                for (int k=0; k<numChildren; k++) {
                    calculate(atomList.getAtom(k), pc);
                }

                PotentialArray intraPotentialArray = getIntraPotentials(targetMolecule.getType());
                final IPotential[] intraPotentials = intraPotentialArray.getPotentials();
                for(int k=0; k<intraPotentials.length; k++) {
                    ((PotentialGroupNbr)intraPotentials[k]).calculateRangeIndependent(targetMolecule,id.direction(), null, pc);
                }
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
    protected void calculate(IAtom atom, PotentialCalculation pc) {
        PotentialArray potentialArray = (PotentialArray)rangedAgentManager.getAgent(atom.getType());
        IPotential[] potentials = potentialArray.getPotentials();
        NeighborCriterion[] criteria = potentialArray.getCriteria();

        for(int i=0; i<potentials.length; i++) {
            switch (potentials[i].nBody()) {
            case 1:
                atomSetSinglet.atom = atom;
                pc.doCalculation(atomSetSinglet, (IPotentialAtomic)potentials[i]);
                break;
            case 2:
                IPotentialAtomic p2 = (IPotentialAtomic)potentials[i];
                NeighborCriterion nbrCriterion = criteria[i];
                neighborIterator.setTarget(atom);
                neighborIterator.reset();
                for (IAtomList pair = neighborIterator.next(); pair != null;
                     pair = neighborIterator.next()) {
                    if (nbrCriterion.accept(pair)) {
                        pc.doCalculation(pair, p2);
                    }
                }
                break;
            }
        }
    }
    
    public static class BoxAgentSiteSource implements BoxAgentSource<BoxCellManager> {
        private final int nCells;
        private final Space space;
        
        public BoxAgentSiteSource(int nCells, Space _space) {
            this.nCells = nCells;
            this.space = _space;
        }

        public BoxCellManager makeAgent(Box box) {
            return new NeighborSiteManager(box,nCells, space);
        }

        public void releaseAgent(BoxCellManager agent) {
        }
    }
}
