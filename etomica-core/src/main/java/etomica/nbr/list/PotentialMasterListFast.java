/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.list;

import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.box.BoxAgentManager;
import etomica.box.BoxCellManager;
import etomica.integrator.Integrator;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculePositionDefinition;
import etomica.nbr.CriterionAdapter;
import etomica.nbr.CriterionSimple;
import etomica.nbr.NeighborCriterion;
import etomica.nbr.PotentialGroupNbr;
import etomica.nbr.cell.NeighborCellManager;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.Debug;

/**
 * PotentialMaster used to implement neighbor listing.  Instance of this
 * class is given as an argument to the Simulation constructor.
 */
public class PotentialMasterListFast extends PotentialMasterList {

    protected final Vector dr;

    /**
     * Constructor specifying space and range for neighbor listing; uses null AtomPositionDefinition.
     */
    public PotentialMasterListFast(Simulation sim, double range, Space _space) {
        this(sim, range, (IMoleculePositionDefinition) null, _space);
    }

    /**
     * Constructs class using given position definition for all atom cell
     * assignments.
     *
     * @param positionDefinition if null, specifies use of atom type's position definition
     */
    public PotentialMasterListFast(Simulation sim, double range, IMoleculePositionDefinition positionDefinition, Space _space) {
        this(sim, range, new BoxAgentSourceCellManagerList(sim, positionDefinition, _space), _space);
    }

    public PotentialMasterListFast(Simulation sim, double range, BoxAgentSourceCellManagerList boxAgentSource, Space _space) {
        this(sim, range, boxAgentSource, new BoxAgentManager<>(boxAgentSource, sim), _space);
    }

    public PotentialMasterListFast(Simulation sim, double range, BoxAgentSourceCellManagerList boxAgentSource, BoxAgentManager<? extends BoxCellManager> agentManager, Space _space) {
        this(sim, range, boxAgentSource, agentManager, new NeighborListAgentSource(range, _space), _space);
    }

    public PotentialMasterListFast(Simulation sim, double range,
                                   BoxAgentSourceCellManagerList boxAgentSource,
                                   BoxAgentManager<? extends BoxCellManager> agentManager,
                                   NeighborListAgentSource neighborListAgentSource, Space _space) {
        super(sim, range, boxAgentSource, agentManager, neighborListAgentSource, _space);
        dr = space.makeVector();
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
        NeighborListManager neighborManager = neighborListAgentManager.getAgent(box);
        if (targetAtom == null && targetMolecule == null) {
            if (Debug.ON && id.direction() != IteratorDirective.Direction.UP) {
                throw new IllegalArgumentException("When there is no target, iterator directive must be up");
            }

            IAtomList atoms = box.getLeafList();
            AtomLeafAgentManager<? extends Integrator.Forcible> agentManager = ((PotentialCalculationForceSum) pc).getAgentManager();
            Boundary boundary = box.getBoundary();
            for (int i = 0; i < atoms.getAtomCount(); i++) {
                IAtom atom = atoms.getAtom(i);
                Vector v = atom.getPosition();
                IAtomList list = neighborManager.getUpList(atom)[0];
                int nNeighbors = list.getAtomCount();
                Integrator.Forcible iAgent = agentManager.getAgent(atom);

                for (int j = 0; j < nNeighbors; j++) {
                    IAtom jAtom = list.getAtom(j);
                    dr.Ev1Mv2(v, jAtom.getPosition());
                    boundary.nearestImage(dr);
                    double r2 = dr.squared();
                    if (r2 > 9) continue;
                    Integrator.Forcible jAgent = agentManager.getAgent(atom);
                    double s6 = 1 / (r2 * r2 * r2);
                    double du = -12 * s6 * (s6 - 0.5);
                    dr.Ea1Tv1(du / r2, dr);
                    jAgent.force().PE(dr);
                    iAgent.force().ME(dr);
                }

            }
        }
        else {

            throw new RuntimeException("not here");
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
    protected void calculate(IMolecule molecule, IteratorDirective.Direction direction, PotentialCalculation pc, NeighborListManager neighborManager) {
        PotentialArray potentialArray = getIntraPotentials(molecule.getType());
        IPotential[] potentials = potentialArray.getPotentials();
        for (int i = 0; i < potentials.length; i++) {
            ((PotentialGroupNbr) potentials[i]).calculateRangeIndependent(molecule, direction, null, pc);
        }

        //cannot use AtomIterator field because of recursive call
        IAtomList list = molecule.getChildList();
        int size = list.getAtomCount();
        for (int i = 0; i < size; i++) {
            calculate(list.getAtom(i), direction, pc, neighborManager);//recursive call
        }
    }

    protected void calculate(IAtom atom, IteratorDirective.Direction direction, PotentialCalculation pc, NeighborListManager neighborManager) {
        singletIterator.setAtom(atom);
        PotentialArray potentialArray = (PotentialArray) rangedAgentManager.getAgent(atom.getType());
        IPotential[] potentials = potentialArray.getPotentials();
        for (int i = 0; i < potentials.length; i++) {
            switch (potentials[i].nBody()) {
                case 1:
                    boolean[] potential1BodyArray = neighborManager.getPotential1BodyList(atom).getInteractingList();
                    if (potential1BodyArray[i]) {
                        atomSetSinglet.atom = atom;
                        pc.doCalculation(atomSetSinglet, (IPotentialAtomic) potentials[i]);
                    }
                    break;
                case 2:
                    if (direction != IteratorDirective.Direction.DOWN) {
                        IAtomList list = neighborManager.getUpList(atom)[i];
                        int nNeighbors = list.getAtomCount();
                        atomPair.atom0 = atom;
                        for (int j = 0; j < nNeighbors; j++) {
                            atomPair.atom1 = list.getAtom(j);
                            pc.doCalculation(atomPair, (IPotentialAtomic) potentials[i]);
                        }
                    }
                    if (direction != IteratorDirective.Direction.UP) {
                        IAtomList list = neighborManager.getDownList(atom)[i];
                        int nNeighbors = list.getAtomCount();
                        atomPair.atom1 = atom;
                        for (int j = 0; j < nNeighbors; j++) {
                            atomPair.atom0 = list.getAtom(j);
                            pc.doCalculation(atomPair, (IPotentialAtomic) potentials[i]);
                        }
                    }
                    break;//switch
                case Integer.MAX_VALUE: //N-body
                    // do the calculation considering the current Atom as the
                    // "central" Atom.
                    if (atomArrayList == null) {
                        atomArrayList = new AtomArrayList();
                    }
                    doNBodyStuff(atom, pc, i, (IPotentialAtomic) potentials[i], neighborManager);
                    if (direction != IteratorDirective.Direction.UP) {
                        // must have a target and be doing "both"
                        // we have to do the calculation considering each of the
                        // target's neighbors
                        IAtomList list = neighborManager.getUpList(atom)[i];
                        for (int j = 0; j < list.getAtomCount(); j++) {
                            IAtom otherAtom = list.getAtom(j);
                            doNBodyStuff(otherAtom, pc, i, (IPotentialAtomic) potentials[i], neighborManager);
                        }
                        list = neighborManager.getDownList(atom)[i];
                        for (int j = 0; j < list.getAtomCount(); j++) {
                            IAtom otherAtom = list.getAtom(j);
                            doNBodyStuff(otherAtom, pc, i, (IPotentialAtomic) potentials[i], neighborManager);
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
        return (NeighborCellManager) boxAgentManager.getAgent(box);
    }

    public int getCellRange() {
        return cellRange;
    }

    public void setCellRange(int newCellRange) {
        cellRange = newCellRange;

        for (BoxCellManager boxCellManager : boxAgentManager.getAgents().values()) {
            ((NeighborCellManager) boxCellManager).setCellRange(cellRange);
        }
    }

}
