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
import etomica.box.BoxCellManager;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.nbr.*;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.species.ISpecies;
import etomica.util.Debug;

import java.util.List;

public class PotentialMasterSite extends PotentialMasterNbr {

    protected final AtomSetSinglet atomSetSinglet;
    protected final AtomsetIteratorPDT neighborIterator;
    private int cellRange;

	/**
	 * Invokes superclass constructor, specifying IteratorFactoryCell
     * for generating molecule iterators.  Sets default nCells of 10 and
     * position definition to null, so that atom type's definition is used
     * to assign cells.
     */
	public PotentialMasterSite(Simulation sim, int nCells, Space space) {
        this(sim, new BoxAgentManager<>(sim, box -> new NeighborSiteManager(box, nCells)), space);
    }

    public PotentialMasterSite(Simulation sim, BoxAgentManager<BoxCellManager> boxAgentManager, Space _space) {
        this(sim, boxAgentManager, new Api1ASite(_space.D(), boxAgentManager));
    }
    
    protected PotentialMasterSite(Simulation sim, BoxAgentManager<? extends BoxCellManager> boxAgentManager, AtomsetIteratorPDT neighborIterator) {
        super(sim, boxAgentManager);
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
        super.addRangedPotentialForTypes(potential, atomType);
        NeighborCriterion criterion;
        if (atomType.length == 2) {
            criterion = new CriterionTypePair(new CriterionAll(), atomType[0], atomType[1]);
            ISpecies moleculeType0 = atomType[0].getSpecies();
            ISpecies moleculeType1 = atomType[1].getSpecies();
            if (moleculeType0 == moleculeType1) {
                criterion = new CriterionInterMolecular(criterion);
            }
            setCriterion(atomType[0], atomType[1], criterion);
        }
        else if (atomType.length == 1) {
            criterion = new CriterionType(new CriterionAll(), atomType[0]);
            setCriterion1Body(potential, atomType[0], criterion);
        }
        else {
            criterion = new CriterionTypesMulti(new CriterionAll(), atomType);
            for (int i = 0; i < atomType.length; i++) {
                for (int j = 0; j <= i; j++) {
                    setCriterion(atomType[i], atomType[j], criterion);
                }
            }
        }
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
        setBoxForCriteria(box);
        IAtom targetAtom = id.getTargetAtom();
        IMolecule targetMolecule = id.getTargetMolecule();
        neighborIterator.setBox(box);
        if (targetAtom == null && targetMolecule == null) {
            if (Debug.ON && id.direction() != IteratorDirective.Direction.UP) {
                throw new IllegalArgumentException("When there is no target, iterator directive must be up");
            }
            neighborIterator.setDirection(IteratorDirective.Direction.UP);
            // invoke setBox on all potentials
            for (IPotential allPotential : allPotentials) {
                allPotential.setBox(box);
            }

            //no target atoms specified
            //call calculate with each molecule
            for (int j=0; j<simulation.getSpeciesCount(); j++) {
                IMoleculeList moleculeList = box.getMoleculeList(simulation.getSpecies(j));
                PotentialArray intraPotentialArray = getIntraPotentials(simulation.getSpecies(j));
                final IPotential[] intraPotentials = intraPotentialArray.getPotentials();
                for (IMolecule molecule : moleculeList) {
                    IAtomList atomList = molecule.getChildList();
                    for (IAtom anAtomList : atomList) {
                        calculate(anAtomList, pc);
                    }

                    for (IPotential intraPotential : intraPotentials) {
                        ((PotentialGroupNbr) intraPotential).calculateRangeIndependent(molecule, id.direction(), null, pc);
                    }
                }
            }
        }
        else {
            // one target atom
            neighborIterator.setDirection(id.direction());
            if (targetAtom != null) {
                IPotentialAtomic[] potentials = getRangedPotentials(targetAtom.getType());
                for (IPotentialAtomic potential : potentials) {
                    if (potential != null) potential.setBox(box);
                }

                //walk up the tree looking for 1-body range-independent potentials that apply to parents
                IMolecule parentMolecule = targetAtom.getParentGroup();
                PotentialArray potentialArray = getIntraPotentials(parentMolecule.getType());
                IPotential[] intraPotentials = potentialArray.getPotentials();
                for (IPotential intraPotential : intraPotentials) {
                    intraPotential.setBox(box);
                    ((PotentialGroupNbr) intraPotential).calculateRangeIndependent(parentMolecule, id.direction(), targetAtom, pc);
                }
                calculate(targetAtom, pc);
            }
            else {
                for (IPotential allPotential : allPotentials) {
                    allPotential.setBox(box);
                }

                IAtomList atomList = (targetMolecule).getChildList();
                for (IAtom anAtomList : atomList) {
                    calculate(anAtomList, pc);
                }

                PotentialArray intraPotentialArray = getIntraPotentials(targetMolecule.getType());
                final IPotential[] intraPotentials = intraPotentialArray.getPotentials();
                for (IPotential intraPotential : intraPotentials) {
                    ((PotentialGroupNbr) intraPotential).calculateRangeIndependent(targetMolecule, id.direction(), null, pc);
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

        List<IPotentialAtomic> potentials1 = rangedPotentials1Body[atom.getType().getIndex()];
        if (potentials1.size() > 0) {
            List<NeighborCriterion> criteria1 = criteria1Body[atom.getType().getIndex()];
            atomSetSinglet.atom = atom;
            for (int i = 0; i < potentials1.size(); i++) {
                if (criteria1.get(i).accept(atom, null)) pc.doCalculation(atomSetSinglet, potentials1.get(i));
            }
        }

        IPotentialAtomic[] potentials = rangedPotentials[atom.getType().getIndex()];
        NeighborCriterion[] myCriteria = criteria[atom.getType().getIndex()];

        for (int i = 0; i < potentials.length; i++) {
            if (potentials[i] == null) continue;
            IPotentialAtomic p2 = potentials[i];
            NeighborCriterion nbrCriterion = myCriteria[i];
            neighborIterator.setTarget(atom);
            neighborIterator.reset();
            for (IAtomList pair = neighborIterator.next(); pair != null;
                 pair = neighborIterator.next()) {
                if (nbrCriterion.accept(pair.get(0), pair.get(1))) {
                    pc.doCalculation(pair, p2);
                }
            }
        }
    }

}
