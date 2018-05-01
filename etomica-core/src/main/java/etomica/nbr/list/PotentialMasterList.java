/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.list;

import etomica.atom.*;
import etomica.box.Box;
import etomica.box.BoxAgentManager;
import etomica.box.BoxCellManager;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.nbr.*;
import etomica.nbr.cell.NeighborCellManager;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.species.ISpecies;
import etomica.util.Debug;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * PotentialMaster used to implement neighbor listing.  Instance of this
 * class is given as an argument to the Simulation constructor.
 */
public class PotentialMasterList extends PotentialMasterNbr {

    private final Space space;
    private final AtomSetSinglet atomSetSinglet;
    private final AtomPair atomPair;
    private final NeighborListAgentSource neighborListAgentSource;
    private final BoxAgentManager<NeighborListManager> neighborListAgentManager;
    private double range;
    private int cellRange;
    private double maxPotentialRange = 0;
    private double safetyFactor = 0.4;

    /**
     * Default constructor uses range of 1.0.
     */
    public PotentialMasterList(Simulation sim, Space _space) {
        this(sim, 1.0, _space);
    }

    public PotentialMasterList(Simulation sim, double range, Space _space) {
        this(sim, range, new NeighborListAgentSource(range), _space);
    }

    public PotentialMasterList(Simulation sim, double range, NeighborListAgentSource neighborListAgentSource, Space _space) {
        super(sim);
        space = _space;
        this.neighborListAgentSource = neighborListAgentSource;
        neighborListAgentSource.setPotentialMaster(this);
        neighborListAgentManager = new BoxAgentManager<>(neighborListAgentSource, sim);
        atomSetSinglet = new AtomSetSinglet();
        atomPair = new AtomPair();
        cellRange = 2;

        // setRange last.  that should always be OK since anyone can call
        // setRange later. if we call it early, member fields won't exist,
        // if we just set range, our stuff associated with existing boxes
        // won't get initialized
        setRange(range);
    }

    /**
     * Convenience method to return the wrapped range-dependent criterion, if
     * one exists
     */
    private static CriterionSimple getRangedCriterion(NeighborCriterion criterion) {
        if (criterion instanceof CriterionSimple) {
            return (CriterionSimple) criterion;
        }
        if (criterion instanceof CriterionAdapter) {
            return getRangedCriterion(((CriterionAdapter) criterion).getWrappedCriterion());
        }
        return null;
    }

    /**
     * Returns the range that determines how far to look for neighbors.
     */
    public double getRange() {
        return range;
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
        recomputeCriteriaRanges();

        neighborListAgentSource.setRange(newRange);

        for (NeighborListManager neighborListManager : neighborListAgentManager.getAgents().values()) {
            neighborListManager.setRange(range);
        }
    }

    /**
     * Returns the maximum range of any potential held by this potential master
     */
    public double getMaxPotentialRange() {
        return maxPotentialRange;
    }

    /**
     * Returns the safety factor.
     */
    public double getSafetyFactor() {
        return safetyFactor;
    }

    /**
     * Sets the safety factor.  The default value, 0.4, is appropriate for most
     * simulations.  Valid values range from 0 to 0.5 (non-inclusive).  The
     * safety factor determines how far an atom can travel from its original
     * position before the neighbor lists are reconstructed.  0.5 means the
     * atom can travel half of its neighbor range.  If another atom also
     * travels half-way then the Atoms could interact without the
     * PotentialMaster recognizing they are neighbors.
     * <p>
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

    public int getCellRange() {
        return cellRange;
    }

    public void setCellRange(int newCellRange) {
        cellRange = newCellRange;

        for (NeighborListManager neighborListManager : neighborListAgentManager.getAgents().values()) {
            neighborListManager.getNeighborCellManager().setCellRange(cellRange);
        }
    }

    public NeighborListManager getNeighborManager(Box box) {
        // we didn't have the simulation when we made the agent manager.
        // setting the simulation after the first time is a quick return
        return neighborListAgentManager.getAgent(box);
    }

    public NeighborCellManager getNbrCellManager(Box box) {
        return neighborListAgentManager.getAgent(box).getNeighborCellManager();
    }

    /**
     * Add the given potential to be used for the given atom types and the
     * given criterion.  If multiple types are given, then the potential will
     * be used for any set of atoms containing only the atom types listed;
     * if types A and B are passed, then the potential will be used for
     * A-A, A-B, B-A and B-B.  The criterion must handle any additional
     * filtering.
     * <p>
     * The given potential will not be held by a PotentialGroup.
     */
    public void addPotentialNbrList(IPotentialAtomic potential, AtomType[] atomTypes, NeighborCriterion criterion) {
        if (potential.getRange() == Double.POSITIVE_INFINITY) {
            throw new RuntimeException("not the method you wanted to call");
        }

        int type1 = atomTypes[0].getIndex();
        int type2 = atomTypes[1].getIndex();
        criteria[type1][type2] = criteria[type2][type1] = criterion;

        if (potential.getRange() > maxPotentialRange) {
            maxPotentialRange = potential.getRange();
        }
        recomputeCriteriaRanges();
    }

    /**
     * Adds the potential as a ranged potential that applies to the given
     * AtomTypes.  This method creates a criterion for the potential and
     * notifies the NeighborListManager of its existence.
     */
    protected void addRangedPotentialForTypes(IPotentialAtomic potential, AtomType[] atomType) {
        super.addRangedPotentialForTypes(potential, atomType);
        // we'll fix the neighbor range later in recomputeCriteriaRanges
        // 0 guarantees the simulation to be hosed if our range is less than the potential range
        // (since recomputeCriteriaRange will bail in that case)
        NeighborCriterion criterion;
        if (potential.nBody() >= 2) {
            NeighborCriterion rangedCriterion;
            if (potential.getRange() < Double.POSITIVE_INFINITY) {
                rangedCriterion = new CriterionSimple(simulation, space, potential.getRange(), 0.0);
            } else {
                // ????? how can this work?
                System.err.println("you have a 'ranged' potential with infinite range.  good luck with that!");
                rangedCriterion = null;
            }
            criterion = rangedCriterion;
            Set<ISpecies> allMySpecies = new HashSet<>();
            allMySpecies.add(atomType[0].getSpecies());
            // if any of our types are in the same species, then enforce inter-molecular
            for (int i = 0; i < atomType.length; i++) {
                if (allMySpecies.contains(atomType[i].getSpecies())) {
                    criterion = new CriterionInterMolecular(criterion);
                    break;
                }
            }
            if (potential.nBody() == 2) {
                setCriterion(atomType[0], atomType[1], criterion);
            } else {
                // n-body, so add criterion for all i-j and also i-i
                for (int i = 0; i < atomType.length; i++) {
                    for (int j = i; j < atomType.length; j++) {
                        setCriterion(atomType[i], atomType[j], criterion);
                    }
                }
            }
        } else {
            setCriterion1Body(potential, atomType[0], new CriterionAll());
        }

        if (potential.getRange() > maxPotentialRange && potential.getRange() < Double.POSITIVE_INFINITY) {
            maxPotentialRange = potential.getRange();
        }
        recomputeCriteriaRanges();
    }

    public void removePotential(IPotentialAtomic potential) {
        super.removePotential(potential);

        maxPotentialRange = 0;
        for (int i = 0; i < rangedPotentials.length; i++) {
            for (int j = 0; j < rangedPotentials.length; j++) {
                if (rangedPotentials[i][j] == null) criteria[i][j] = null;
                double pRange = rangedPotentials[i][j].getRange();
                if (pRange == Double.POSITIVE_INFINITY) {
                    continue;
                }
                if (pRange > maxPotentialRange) {
                    maxPotentialRange = pRange;
                }
            }
        }

        recomputeCriteriaRanges();
    }

    @Override
    public BoxCellManager getBoxCellManager(Box box) {
        return neighborListAgentManager.getAgent(box).getNeighborCellManager();
    }

    /**
     * Recomputes the maximum potential range (which might change without this
     * class receiving notification) and readjust cell lists
     */
    public void reset() {
        maxPotentialRange = 0;
        for (int i = 0; i < rangedPotentials.length; i++) {
            for (int j = i; j < rangedPotentials.length; j++) {
                IPotentialAtomic p = rangedPotentials[i][j];
                if (p != null && p.getRange() > maxPotentialRange) {
                    maxPotentialRange = p.getRange();
                }
            }
        }
        recomputeCriteriaRanges();

        neighborListAgentManager.getAgents().values().forEach(NeighborListManager::reset);
    }

    /**
     * Recomputes the range for all criterion based on our own range.  The
     * range for each criterion is set so that they all have an equal max
     * displacement.  Our nominal neighbor range is used for the criterion
     * with the longest potential range.
     */
    private void recomputeCriteriaRanges() {
        double maxDisplacement = (getRange() - maxPotentialRange) * safetyFactor;
        if (maxDisplacement < 0) {
            // someone probably added a long ranged potential and hasn't updated the PotentialMaster's range
            // when they do, we'll get called again.
            // if they don't, the simulation will probably crash
            return;
        }
        for (int i = 0; i < rangedPotentials.length; i++) {
            for (int j = i; j < rangedPotentials.length; j++) {
                IPotentialAtomic p = rangedPotentials[i][j];
                NeighborCriterion c = criteria[i][j];
                if (c == null) continue;
                CriterionSimple rangedCriterion = getRangedCriterion(c);
                if (rangedCriterion != null) {
                    double newRange = maxDisplacement / safetyFactor + p.getRange();
                    rangedCriterion.setNeighborRange(newRange);
                    rangedCriterion.setInteractionRange(p.getRange());
                    rangedCriterion.setSafetyFactor(safetyFactor);
                }
            }
        }
    }

    public void calculate(Box box, PotentialCalculation pc, boolean includeLrc) {
        // invoke setBox on all potentials
        setBoxForPotentials(box);
        NeighborListManager nbrManager = neighborListAgentManager.getAgent(box);
        IAtomList atoms = box.getLeafList();
        for (int i = 0; i < atoms.size(); i++) {
            calculate(atoms.get(i), IteratorDirective.Direction.UP, pc, nbrManager);
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
        if (!enabled) return;
        IAtom targetAtom = id.getTargetAtom();
        IMolecule targetMolecule = id.getTargetMolecule();
        if (targetAtom == null && targetMolecule == null) {
            if (Debug.ON && id.direction() != IteratorDirective.Direction.UP) {
                throw new IllegalArgumentException("When there is no target, iterator directive must be up");
            }

            calculate(box, pc, id.includeLrc);
        } else {
            NeighborListManager neighborManager = neighborListAgentManager.getAgent(box);
            if (id.direction() != IteratorDirective.Direction.UP) {
                neighborManager.ensureDownLists();
            }

            if (targetAtom != null) {
                int targetTypeIdx = targetAtom.getType().getIndex();
                for (int j = 0; j < rangedPotentials[targetTypeIdx].length; j++) {
                    if (rangedPotentials[targetTypeIdx][j] != null) rangedPotentials[targetTypeIdx][j].setBox(box);
                }

                //first walk up the tree looking for 1-body range-independent potentials that apply to parents
                IMolecule parentAtom = targetAtom.getParentGroup();
                PotentialArray potentialArray = getIntraPotentials(parentAtom.getType());
                IPotential[] potentials = potentialArray.getPotentials();
                for (int i = 0; i < potentials.length; i++) {
                    potentials[i].setBox(box);
                    ((PotentialGroupNbr) potentials[i]).calculateRangeIndependent(parentAtom, id.direction(), targetAtom, pc);
                }
                calculate(targetAtom, id.direction(), pc, neighborManager);
            } else {
                PotentialArray potentialArray = getIntraPotentials(targetMolecule.getType());
                IPotential[] potentials = potentialArray.getPotentials();
                for (int i = 0; i < potentials.length; i++) {
                    potentials[i].setBox(box);
                }

                calculate(targetMolecule, id.direction(), pc, neighborManager);
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
    private void calculate(IMolecule molecule, IteratorDirective.Direction direction, PotentialCalculation pc, NeighborListManager neighborManager) {
        PotentialArray potentialArray = getIntraPotentials(molecule.getType());
        IPotential[] potentials = potentialArray.getPotentials();
        for (int i = 0; i < potentials.length; i++) {
            ((PotentialGroupNbr) potentials[i]).calculateRangeIndependent(molecule, direction, null, pc);
        }

        //cannot use AtomIterator field because of recursive call
        IAtomList list = molecule.getChildList();
        int size = list.size();
        for (int i = 0; i < size; i++) {
            calculate(list.get(i), direction, pc, neighborManager);//recursive call
        }
    }

    private void calculate(IAtom atom, IteratorDirective.Direction direction, PotentialCalculation pc, NeighborListManager neighborManager) {
        List<IPotentialAtomic> potentials1 = rangedPotentials1Body[atom.getType().getIndex()];
        if (!potentials1.isEmpty()) {
            boolean[] potential1BodyArray = neighborManager.getPotential1BodyList(atom).getInteractingList();
            atomSetSinglet.atom = atom;
            for (int i = 0; i < potentials1.size(); i++) {
                if (potential1BodyArray[i]) {
                    atomSetSinglet.atom = atom;
                    pc.doCalculation(atomSetSinglet, potentials1.get(i));
                }
            }
        }

        IPotentialAtomic[] potentials = rangedPotentials[atom.getType().getIndex()];
        for (int i = 0; i < potentials.length; i++) {
            if (potentials[i] == null) continue;
            switch (potentials[i].nBody()) {
                case 2:
                    if (direction != IteratorDirective.Direction.DOWN) {
                        IAtomList list = neighborManager.getUpList(atom)[i];
                        int nNeighbors = list.size();
                        atomPair.atom0 = atom;
                        for (int j = 0; j < nNeighbors; j++) {
                            atomPair.atom1 = list.get(j);
                            pc.doCalculation(atomPair, potentials[i]);
                        }
                    }
                    if (direction != IteratorDirective.Direction.UP) {
                        IAtomList list = neighborManager.getDownList(atom)[i];
                        int nNeighbors = list.size();
                        atomPair.atom1 = atom;
                        for (int j = 0; j < nNeighbors; j++) {
                            atomPair.atom0 = list.get(j);
                            pc.doCalculation(atomPair, potentials[i]);
                        }
                    }
                    break;//switch
                case Integer.MAX_VALUE: //N-body
                    calculateNBody(atom, direction, pc, neighborManager, potentials[i], i);


            }
        }
    }

    private static void calculateNBody(IAtom atom, IteratorDirective.Direction direction, PotentialCalculation pc, NeighborListManager neighborManager, IPotentialAtomic potential, int i) {
        // do the calculation considering the current Atom as the
        // "central" Atom.
        doNBodyStuff(atom, pc, i, potential, neighborManager);
        if (direction != IteratorDirective.Direction.UP) {
            // must have a target and be doing "both"
            // we have to do the calculation considering each of the
            // target's neighbors
            IAtomList upList = neighborManager.getUpList(atom)[i];
            for (int j = 0; j < upList.size(); j++) {
                IAtom otherAtom = upList.get(j);
                doNBodyStuff(otherAtom, pc, i, potential, neighborManager);
            }
            IAtomList downList = neighborManager.getDownList(atom)[i];
            for (int j = 0; j < downList.size(); j++) {
                IAtom otherAtom = downList.get(j);
                doNBodyStuff(otherAtom, pc, i, potential, neighborManager);
            }
        }
    }

    /**
     * Invokes the PotentialCalculation for the given Atom with its up and down
     * neighbors as a single AtomSet.
     */
    private static void doNBodyStuff(IAtom atom, PotentialCalculation pc, int potentialIndex,
                              IPotentialAtomic potential, NeighborListManager neighborManager) {
        IAtomList upList = neighborManager.getUpList(atom)[potentialIndex];
        IAtomList downList = neighborManager.getDownList(atom)[potentialIndex];

        // Pre-sizing avoids array reallocs in addAll, it can simply do System.arrayCopy each time
        AtomArrayList atomList = new AtomArrayList(upList.size() + downList.size() + 1);
        atomList.add(atom);
        atomList.addAll(upList);
        atomList.addAll(downList);
        pc.doCalculation(atomList, potential);
    }

}
