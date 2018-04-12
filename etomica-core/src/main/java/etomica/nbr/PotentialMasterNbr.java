/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr;

import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.box.BoxAgentManager;
import etomica.box.BoxCellManager;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.species.ISpecies;

import java.util.ArrayList;
import java.util.List;

public abstract class PotentialMasterNbr extends PotentialMaster {

    protected final IPotentialAtomic[][] rangedPotentials;
    protected final List<IPotentialAtomic>[] rangedPotentials1Body;
    protected final List<IPotentialAtomic> potentials0body;
    protected final NeighborCriterion[][] criteria;
    protected final List<NeighborCriterion>[] criteria1Body;
    private final PotentialArray[] intraPotentials;
    protected final Simulation simulation;

    public PotentialMasterNbr(Simulation sim) {
        super();
        simulation = sim;
        int numAtomTypes = 0;
        if (sim.getSpeciesCount() > 0) {
            ISpecies lastSpecies = sim.getSpecies(sim.getSpeciesCount() - 1);
            numAtomTypes = lastSpecies.getAtomType(lastSpecies.getAtomTypeCount() - 1).getIndex() + 1;
        }
        rangedPotentials = new IPotentialAtomic[numAtomTypes][numAtomTypes];
        criteria = new NeighborCriterion[numAtomTypes][numAtomTypes];

        rangedPotentials1Body = new ArrayList[numAtomTypes];
        criteria1Body = new ArrayList[numAtomTypes];
        for (int i = 0; i < numAtomTypes; i++) {
            rangedPotentials1Body[i] = new ArrayList<>();
            criteria1Body[i] = new ArrayList<>();
        }

        potentials0body = new ArrayList<IPotentialAtomic>();

        intraPotentials = new PotentialArray[sim.getSpeciesList().size()];
        for (int i = 0; i < intraPotentials.length; i++) {
            intraPotentials[i] = new PotentialArray();
        }
    }

    public PotentialGroup makePotentialGroup(int nBody) {
        return new PotentialGroupNbr(nBody);
    }

    public void addPotential(IPotentialMolecular potential, ISpecies[] species) {
        if (!(potential instanceof PotentialGroup)) {
            System.err.println("You gave me a concrete molecule potential and I'm very confused now.  I'll pretend like that's OK but don't hold your breath.");
        }
        super.addPotential(potential, species);
    }

    public void potentialAddedNotify(IPotentialAtomic subPotential, PotentialGroup pGroup) {
        super.potentialAddedNotify(subPotential, pGroup);
        AtomType[] atomTypes = pGroup.getAtomTypes(subPotential);
        if (atomTypes == null) {
            //change
            if (pGroup.nBody() < 2 && subPotential.getRange() == Double.POSITIVE_INFINITY) {
                //pGroup is PotentialGroupNbr
                if (subPotential.nBody() == 0) {
                    potentials0body.add(subPotential);
                } else if (pGroup.nBody() == 1) {
                    ISpecies[] parentType = getSpecies(pGroup);
                    intraPotentials[parentType[0].getIndex()].addPotential(pGroup);
                }
            } else {
                //FIXME what to do with this case?  Fail!
                System.err.println("You have a child-potential of a 2-body PotentialGroup or range-dependent potential, but it's not type-based.  Enjoy crashing or fix bug 85");
            }
            return;
        }
        addRangedPotentialForTypes(subPotential, atomTypes);
    }

    protected void addRangedPotentialForTypes(IPotentialAtomic subPotential, AtomType[] atomTypes) {
        if (atomTypes.length == 1) {
            if (subPotential.nBody() == 1) {
                rangedPotentials1Body[atomTypes[0].getIndex()].add(subPotential);
            } else {
                rangedPotentials[atomTypes[0].getIndex()][atomTypes[0].getIndex()] = subPotential;
            }
            return;
        }
        for (int i = 0; i < atomTypes.length; i++) {
            int ti = atomTypes[i].getIndex();
            for (int j = i + 1; j < atomTypes.length; j++) {
                int tj = atomTypes[j].getIndex();
                rangedPotentials[ti][tj] = rangedPotentials[tj][ti] = subPotential;
            }
        }
    }

    public void removePotential(IPotentialAtomic potential) {
        super.removePotential(potential);
        if (potential.getRange() < Double.POSITIVE_INFINITY) {
            if (potential.nBody() == 1) {
                for (int i = 0; i < rangedPotentials1Body.length; i++) {
                    rangedPotentials1Body[i].remove(potential);
                }
            } else {
                for (int i = 0; i < rangedPotentials.length; i++) {
                    for (int j = 0; j < rangedPotentials.length; j++) {
                        if (rangedPotentials[i][j] == potential) rangedPotentials[i][j] = null;
                    }
                }
            }
        } else if (potential instanceof PotentialGroup) {
            for (PotentialArray potentialArray : this.intraPotentials) {
                potentialArray.removePotential(potential);
            }
        }
    }

    public final IPotentialAtomic[][] getRangedPotentials() {
        return rangedPotentials;
    }

    public final IPotentialAtomic[] getRangedPotentials(AtomType atomType) {
        return rangedPotentials[atomType.getIndex()];
    }

    public final List<IPotentialAtomic> getRangedPotentials1Body(AtomType atomType) {
        return rangedPotentials1Body[atomType.getIndex()];
    }

    public final PotentialArray getIntraPotentials(ISpecies species) {
        return intraPotentials[species.getIndex()];
    }


    /**
     * @param type1 the first atom type
     * @param type2 the second atom type
     * @return the criterion used to determine when atoms of the given type interact.
     */
    public NeighborCriterion getCriterion(AtomType type1, AtomType type2) {
        return criteria[type1.getIndex()][type2.getIndex()];
    }

    /**
     * @param type the atom type for which you want to know the NeighborCriteria
     * @return an array of NeighborCrtieria for the given AtomType with each
     * other AtomType.
     */
    public NeighborCriterion[] getCriteria(AtomType type) {
        return criteria[type.getIndex()];
    }

    /**
     * @param type the atom type for which you want to know the NeighborCriteria
     * @return aa List of 1-body NeighborCrtieria for the given AtomType with each
     * other AtomType.
     */
    public List<NeighborCriterion> getCriteria1Body(AtomType type) {
        return criteria1Body[type.getIndex()];
    }

    /**
     * Sets the criterion associated with the given pair of atom types, overriding
     * the default provided by the PotentialMasterList.  The criterion can be
     * configured by calling getCriterion(Potential) and changing the
     * criterion.
     */
    public void setCriterion(AtomType type1, AtomType type2, NeighborCriterion criterion) {
        int t1 = type1.getIndex();
        int t2 = type2.getIndex();
        criteria[t1][t2] = criteria[t2][t1] = criterion;
    }

    /**
     * Sets the criterion associated with the given 1-body potential, overriding
     * the default provided by the PotentialMasterList.  The criterion can be
     * configured by calling getCriterion(Potential) and changing the
     * criterion.  The potential passed to this method must be a potential
     * handled by this instance.
     */
    public void setCriterion1Body(IPotentialAtomic potential, AtomType type, NeighborCriterion criterion) {
        List<IPotentialAtomic> potentials = rangedPotentials1Body[type.getIndex()];
        List<NeighborCriterion> myCriteria = criteria1Body[type.getIndex()];
        int idx = potentials.indexOf(potential);
        if (idx < 0) throw new RuntimeException("Could not find potential in my list");
        // remove the criterion to all existing NeighborListManagers
        if (idx < myCriteria.size()) {
            myCriteria.set(idx, criterion);
        } else {
            myCriteria.add(criterion);
        }
    }

    /**
     * Sets the box for all NeighborCriterions
     *
     * @param box the box!
     */
    public void setBoxForCriteria(Box box) {
        for (int i = 0; i < criteria.length; i++) {
            for (int j = 0; j <= i; j++) {
                if (criteria[i][j] != null) criteria[i][j].setBox(box);
            }
        }
        for (int i = 0; i < criteria1Body.length; i++) {
            for (NeighborCriterion c : criteria1Body[i]) {
                c.setBox(box);
            }
        }
    }

    /**
     * Sets the box for all ranged potentials
     *
     * @param box the box!
     */
    protected void setBoxForPotentials(Box box) {
        for (int i = 0; i < rangedPotentials.length; i++) {
            for (int j = 0; j <= i; j++) {
                if (rangedPotentials[i][j] != null) rangedPotentials[i][j].setBox(box);
            }
        }
        for (int i = 0; i < rangedPotentials1Body.length; i++) {
            for (IPotentialAtomic p : rangedPotentials1Body[i]) {
                p.setBox(box);
            }
        }
        for (IPotentialAtomic p : potentials0body) {
            p.setBox(box);
        }
        for (PotentialArray pa : intraPotentials) {
            IPotential[] allp = pa.getPotentials();
            for (IPotential p : allp) {
                p.setBox(box);
            }
        }
    }

    public abstract BoxCellManager getBoxCellManager(Box box);
}
