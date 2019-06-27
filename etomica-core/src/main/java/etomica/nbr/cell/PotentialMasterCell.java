/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.cell;

import etomica.atom.*;
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

/**
 * A PotentialMaster for use with a Simulation where cell-listing of atoms and
 * their neighbors is appropriate.
 *
 * @author Andrew Schultz
 */
public class PotentialMasterCell extends PotentialMasterNbr {

    private double range;
    private int cellRange;
    private final BoxAgentSourceCellManager cellManagerSource;
    private final BoxAgentManager<NeighborCellManager> neighborCellManagers;
    private final BoxAgentManager<NeighborIterator> neighborIterators;

    /**
     * Creates PotentialMasterCell with default (1.0) range.  Range
     * should be set manually via setRange method.
     */
    public PotentialMasterCell(Simulation sim, Space _space) {
        this(sim, 1.0, _space);
    }

    /**
     * Constructs with null AtomPositionDefinition, which indicates the position
     * definition given with each atom's AtomType should be used.
     *
     * @param _space the governing Space
     * @param range  the neighbor distance.  May be changed after construction.
     */
    public PotentialMasterCell(Simulation sim, double range, Space _space) {
        this(sim, range, new BoxAgentSourceCellManager(null, range));
    }

    public PotentialMasterCell(Simulation sim, double range, BoxAgentSourceCellManager boxAgentSource) {
        super(sim);
        this.range = range;
        this.cellManagerSource = boxAgentSource;
        this.neighborCellManagers = new BoxAgentManager<>(boxAgentSource, sim);
        this.neighborIterators = new BoxAgentManager<NeighborIterator>(sim, box -> new NeighborIteratorCell(neighborCellManagers.getAgent(box)));
        setRange(range);
    }

    public double getRange() {
        return range;
    }

    public void setRange(double d) {
        cellManagerSource.setRange(d);
        range = d;

        for (NeighborCellManager manager : neighborCellManagers.getAgents().values()) {
            manager.setPotentialRange(range);
        }
    }

    public void setCellRange(int d) {
        this.cellRange = d;

        for (NeighborCellManager neighborCellManager : this.neighborCellManagers.getAgents().values()) {
            neighborCellManager.setCellRange(this.cellRange);
        }
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
        return this.neighborCellManagers.getAgent(box);
    }

    public NeighborCellManager getNbrCellManager(Box box) {
        return this.neighborCellManagers.getAgent(box);

        // TODO: why??
//        NeighborCellManager manager = boxAgentManagerNeighborCell.getAgent(box);
//        manager.setPotentialRange(range);
//        int cr = getCellRange();
//        if (cr == 0) throw new RuntimeException("need to set cell range first");
//        manager.setCellRange(cr);
//        return manager;
    }

    /**
     * Reassign atoms to cell lists for all boxes.
     */
    public void reset() {
        this.neighborCellManagers.getAgents().values().forEach(BoxCellManager::assignCellAll);
    }

    public void calculate(Box box, PotentialCalculation pc, boolean includeLrc) {

        setBoxForCriteria(box);
        setBoxForPotentials(box);

        IAtomList atoms = box.getLeafList();
        for (int i = 0; i < atoms.size(); i++) {
            calculate(atoms.get(i), neighborIterators.getAgent(box), pc, IteratorDirective.Direction.UP);
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
                setBoxForCriteria(box);
                setBoxForPotentials(box);

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
        NeighborCriterion[] myCriteria = criteria[atom.getType().getIndex()];
        IPotentialAtomic[] potentials = rangedPotentials[atom.getType().getIndex()];

        neighborIterator.forEachNeighbor(atom, direction,
                (targetAtom, otherAtom) -> {
                    NeighborCriterion criterion = myCriteria[otherAtom.getType().getIndex()];
                    if (criterion != null && criterion.accept(targetAtom, otherAtom)) {
                        pc.doCalculation(new AtomPair(targetAtom, otherAtom), potentials[otherAtom.getType().getIndex()]);
                    }
                },
                (otherAtom, targetAtom) -> {
                    NeighborCriterion criterion = myCriteria[otherAtom.getType().getIndex()];
                    if (criterion != null && criterion.accept(otherAtom, targetAtom)) {
                        pc.doCalculation(new AtomPair(otherAtom, targetAtom), potentials[otherAtom.getType().getIndex()]);
                    }
                }
        );

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
