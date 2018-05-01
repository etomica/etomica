/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.site;

import etomica.atom.*;
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
import java.util.function.Function;

public class PotentialMasterSite extends PotentialMasterNbr {

    private final BoxAgentManager<NeighborIterator> neighborIterators;
    private int cellRange;
    private final BoxAgentManager<NeighborSiteManager> neighborSiteManagers;

	/**
	 * Invokes superclass constructor, specifying IteratorFactoryCell
     * for generating molecule iterators.  Sets default nCells of 10 and
     * position definition to null, so that atom type's definition is used
     * to assign cells.
     */
	public PotentialMasterSite(Simulation sim, int nCells, Space space) {
        this(sim, nCells);
    }

    protected PotentialMasterSite(Simulation sim, int nCells) {
        super(sim);
        this.neighborSiteManagers = new BoxAgentManager<>(sim, box -> new NeighborSiteManager(box, nCells));
        this.neighborIterators = new BoxAgentManager<>(sim, box -> new NeighborIteratorSite(neighborSiteManagers.getAgent(box), box));
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
            criterion = new CriterionAll();
            ISpecies moleculeType0 = atomType[0].getSpecies();
            ISpecies moleculeType1 = atomType[1].getSpecies();
            if (moleculeType0 == moleculeType1) {
                criterion = new CriterionInterMolecular(criterion);
            }
            setCriterion(atomType[0], atomType[1], criterion);
        }
        else if (atomType.length == 1) {
            setCriterion1Body(potential, atomType[0], new CriterionAll());
        }
        else {
            criterion = new CriterionAll();
            for (int i = 0; i < atomType.length; i++) {
                for (int j = 0; j <= i; j++) {
                    setCriterion(atomType[i], atomType[j], criterion);
                }
            }
        }
    }

    @Override
    public BoxCellManager getBoxCellManager(Box box) {
        return neighborSiteManagers.getAgent(box);
    }


    public void calculate(Box box, PotentialCalculation pc, boolean includeLrc) {

        NeighborIterator nbrIterator = neighborIterators.getAgent(box);
        setBoxForCriteria(box);
        setBoxForPotentials(box);

        IAtomList atoms = box.getLeafList();
        for (int i = 0; i < atoms.size(); i++) {
            calculate(atoms.get(i), nbrIterator, pc, IteratorDirective.Direction.UP);
        }
        for (int i = 0; i < simulation.getSpeciesCount(); i++) {
            PotentialArray intraPotentialArray = getIntraPotentials(simulation.getSpecies(i));
            IPotential[] intraPotentials = intraPotentialArray.getPotentials();
            if (intraPotentials.length > 0) {
                IMoleculeList moleculeList = box.getMoleculeList(simulation.getSpecies(i));
                for (int j = 0; j < moleculeList.size(); j++) {
                    for (IPotential intraPotential : intraPotentials) {
                        ((PotentialGroupNbr) intraPotential).calculateRangeIndependent(moleculeList.get(i), IteratorDirective.Direction.UP, null, pc);
                    }
                }
            }
        }
    }

    @Override
    public void calculate(Box box, IteratorDirective id, PotentialCalculation pc) {
        if (!enabled) {
            return;
        }

        IAtom targetAtom = id.getTargetAtom();
        IMolecule targetMolecule = id.getTargetMolecule();
        if (targetAtom == null && targetMolecule == null) {
            if (Debug.ON && id.direction() != IteratorDirective.Direction.UP) {
                throw new IllegalArgumentException("When there is no target, iterator directive must be up");
            }

            calculate(box, pc, id.includeLrc);

        } else {
            if (targetAtom != null) {
                for (IPotentialAtomic potential : getRangedPotentials(targetAtom.getType())) {
                    if (potential != null) {
                        potential.setBox(box);
                    }
                }


                //look for 1-body range-independent potentials that apply to parents
                IMolecule parentMolecule = targetAtom.getParentGroup();
                PotentialArray potentialArray = getIntraPotentials(parentMolecule.getType());
                IPotential[] intraPotentials = potentialArray.getPotentials();
                for (IPotential intraPotential : intraPotentials) {
                    intraPotential.setBox(box);
                    ((PotentialGroupNbr) intraPotential).calculateRangeIndependent(parentMolecule, id.direction(), targetAtom, pc);
                }
                calculate(targetAtom, neighborIterators.getAgent(box), pc, id.direction());
            } else {
                for (int i = 0; i < targetMolecule.getChildList().size(); i++) {
                    calculate(targetMolecule.getChildList().get(i), neighborIterators.getAgent(box), pc, id.direction());
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

    private void calculate(IAtom atom, NeighborIterator neighborIterator, PotentialCalculation pc, IteratorDirective.Direction direction) {
        calculate1Body(atom, pc);

        IPotentialAtomic[] potentials = rangedPotentials[atom.getType().getIndex()];
        NeighborCriterion[] myCriteria = criteria[atom.getType().getIndex()];

        neighborIterator.forEachNeighbor(atom, direction, (atom1, atom2) -> {
            for (int i = 0; i < potentials.length; i++) {
                if (potentials[i] == null) {
                    continue;
                }

                IPotentialAtomic p2 = potentials[i];
                NeighborCriterion nbrCriterion = myCriteria[i];
                if (nbrCriterion.accept(atom1, atom2)) {
                    pc.doCalculation(new AtomPair(atom1, atom2), p2);
                }
            }
        });

    }

    private void calculate1Body(IAtom atom, PotentialCalculation pc) {
        List<IPotentialAtomic> potentials1 = rangedPotentials1Body[atom.getType().getIndex()];
        if (!potentials1.isEmpty()) {
            List<NeighborCriterion> criteria1 = criteria1Body[atom.getType().getIndex()];
            for (int i = 0; i < potentials1.size(); i++) {
                if (criteria1.get(i).accept(atom, null)) {
                    pc.doCalculation(new AtomSetSinglet(atom), potentials1.get(i));
                }
            }
        }
    }

}
