/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.api.ISpecies;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.iterator.AtomsetIteratorBasisDependent;
import etomica.box.Box;
import etomica.chem.models.Model;
import etomica.chem.models.Model.PotentialAndIterator;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeArrayList;
import etomica.molecule.iterator.IteratorFactory;
import etomica.molecule.iterator.MoleculeIteratorAll;
import etomica.molecule.iterator.MoleculesetIteratorPDT;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


/**
 * Manager of all potentials in simulation.
 * Most calls to compute the energy or other potential calculations begin
 * with the calculate method of this class.  It then passes the calculation
 * on to the contained potentials.
 *
 * @author David Kofke
 */
public class PotentialMaster {

    protected PotentialMasterLrc lrcMaster;
    protected IteratorFactory iteratorFactory;
    protected List<PotentialLinker> potentialList;
    protected boolean enabled = true;

    public PotentialMaster() {
        this(IteratorFactory.INSTANCE);
    }

    public PotentialMaster(IteratorFactory iteratorFactory) {
        this.iteratorFactory = iteratorFactory;
        potentialList = new ArrayList<PotentialLinker>();
    }

    /**
     * Returns the object that oversees the long-range
     * correction zero-body potentials.
     */
    public PotentialMasterLrc lrcMaster() {
        if (lrcMaster == null) lrcMaster = new PotentialMasterLrc();
        return lrcMaster;
    }

    /**
     * Returns an nBody PotentialGroup appropriate for this type of
     * PotentialMaster.
     */
    public PotentialGroup makePotentialGroup(int nBody) {
        return new PotentialGroup(nBody);
    }

    /**
     * Performs the given PotentialCalculation on the atoms of the given Box.
     * Sets the box for all molecule iterators and potentials, sets target
     * and direction for iterators as specified by given IteratorDirective,
     * and applies doCalculation of given PotentialCalculation with the iterators
     * and potentials.
     */
    public void calculate(Box box, IteratorDirective id, PotentialCalculation pc) {
        if (!enabled) return;
        IMolecule targetMolecule = id.getTargetMolecule();
        IAtom targetAtomLeaf = id.getTargetAtom();
        if (targetAtomLeaf != null) {
            targetMolecule = targetAtomLeaf.getParentGroup();
        }

        for (PotentialLinker link : potentialList) {
            if (!link.enabled) continue;
            final MoleculesetIteratorPDT atomIterator = link.iterator;
            final IPotentialMolecular potential = link.potential;
            atomIterator.setBox(box);
            potential.setBox(box);
            atomIterator.setTarget(targetMolecule);
            atomIterator.setDirection(id.direction());
            if (pc instanceof PotentialCalculationMolecular) {
                atomIterator.reset();
                for (IMoleculeList atoms = atomIterator.next(); atoms != null;
                     atoms = atomIterator.next()) {
                    ((PotentialCalculationMolecular) pc).doCalculation(atoms, potential);
                }
            } else if (potential instanceof PotentialGroup) {
                ((PotentialGroup) potential).calculate(atomIterator, id.direction(), targetAtomLeaf, pc);
            }
        }

        if (lrcMaster != null) {
            lrcMaster.calculate(box, id, pc);
        }
    }

    /**
     * Add the given Model's intramolecular potentials to this PotentialMaster
     */
    public void addModel(Model newModel) {
        if (getPotential(new ISpecies[]{newModel.getSpecies()}) != null) {
            throw new IllegalArgumentException(newModel + " has already been added");
        }
        PotentialAndIterator[] potentialsAndIterators = newModel.getPotentials();
        PotentialGroup pGroup = makePotentialGroup(1);
        for (int i = 0; i < potentialsAndIterators.length; i++) {
            pGroup.addPotential(potentialsAndIterators[i].getPotential(),
                    potentialsAndIterators[i].getIterator());
        }
        addPotential(pGroup, new ISpecies[]{newModel.getSpecies()});
    }

    /**
     * Indicates to the PotentialMaster that the given potential should apply to
     * the specified species.  Exception is thrown if the potential.nBody() value
     * is different from the length of the species array.  Thus, for example, if
     * giving a 2-body potential, then the array should contain exactly
     * two species; the species may refer to the same instance (appropriate for an
     * intra-species potential, defining the interactions between molecules of the
     * same species).
     */
    public void addPotential(IPotentialMolecular potential, ISpecies[] species) {
        if (potential.nBody() == 0) {
            addPotential(potential, new MoleculeIterator0(), null);
        } else if (potential.nBody() == Integer.MAX_VALUE) {
            addPotential(potential, new MoleculeIteratorAll(species), null);
        } else if (potential.nBody() != species.length) {
            throw new IllegalArgumentException("Illegal species length");
        } else {
            MoleculesetIteratorPDT iterator = iteratorFactory.makeMoleculeIterator(species);
            addPotential(potential, iterator, species);
        }
    }

    /**
     * Indicates to the PotentialMaster that the given potential should apply to
     * the specified atom types.  The potential is assumed to be intermolecular.
     * The given types should not include any type which is the descendent of
     * another.  Potential group hierarchy will be constructed as needed above
     * the level of the given atom types.
     */
    public void addPotential(IPotentialAtomic potential, AtomType[] atomTypes) {
        if (potential.nBody() != Integer.MAX_VALUE && potential.nBody() != atomTypes.length) {
            throw new IllegalArgumentException("nBody of potential must match number of atom types");
        }

        // depth of molecules
        ISpecies[] parentAtomTypes = new ISpecies[atomTypes.length];
        for (int i = 0; i < atomTypes.length; i++) {
            parentAtomTypes[i] = atomTypes[i].getSpecies();
        }
        // look for a PotentialGroup that applies to parentAtomTypes
        PotentialGroup pGroup = getPotential(parentAtomTypes);
        if (pGroup == null) { // didn't find an appropriate potentialgroup
            pGroup = makePotentialGroup(potential.nBody());
            addPotential(pGroup, parentAtomTypes);
        }
        pGroup.addPotential(potential, atomTypes);
    }

    /**
     * Notifies the PotentialMaster that the sub-potential has been added to
     * the given PotentialGroup, which is associated (but not necessarily held
     * by) this PotentialMaster.
     * This method is called by PotentialGroup and should not be called in
     * other circumstances.
     */
    public void potentialAddedNotify(IPotentialAtomic subPotential, PotentialGroup pGroup) {
        // do nothing.  this is here for subclasses to override
    }

    /**
     * Returns the potential that applies to the specified types,
     * or null of no existing potential applies.
     */
    public PotentialGroup getPotential(ISpecies[] types) {
        for (PotentialLinker link : potentialList) {
            if (link.potential instanceof PotentialGroup) {
                if (Arrays.equals(types, link.types)) {
                    return (PotentialGroup) link.potential;
                }
            }
        }
        return null;
    }

    public ISpecies[] getSpecies(IPotentialMolecular potential) {
        for (PotentialLinker link : potentialList) {
            if (link.potential == potential) {
                return link.types;
            }
        }
        return null;
    }

    public void addPotential(IPotentialMolecular potential, MoleculesetIteratorPDT iterator, ISpecies[] types) {
        //the order of the given potential should be consistent with the order of the iterator
        if (potential.nBody() != iterator.nBody()) {
            throw new RuntimeException("Error: adding to PotentialGroup a potential and iterator that are incompatible");
        }
        //Set up to evaluate zero-body potentials last, since they may need other potentials
        //to be configured for calculation first
        potentialList.add(new PotentialLinker(potential, iterator, types));
        if (potential instanceof PotentialGroup) {
            ((PotentialGroup) potential).setPotentialMaster(this);
        }
    }

    /**
     * Removes given potential from the group.  No error is generated if
     * potential is not in group.
     */
    public synchronized void removePotential(IPotentialMolecular potential) {

        for (PotentialLinker link : potentialList) {
            if (link.potential == potential) {
                potentialList.remove(link);
                return;
            }
        }
        throw new RuntimeException("potential not found");
    }

    /**
     * Removes given potential from the group.  No error is generated if
     * potential is not in group.
     */
    public synchronized void removePotential(IPotentialAtomic potential) {

        for (PotentialLinker link : potentialList) {
            if (link.potential instanceof PotentialGroup && ((PotentialGroup) link.potential).contains(potential)) {
                ((PotentialGroup) link.potential).removePotential(potential);
                return;
            }
        }
        throw new RuntimeException("potential not found");
    }

    /**
     * @return Returns enabled flag.
     */
    public boolean isEnabled() {
        return enabled;
    }

    /**
     * Permits enabling/disabling of all potentials.  Default is enabled (true).
     *
     * @param enabled flags if potentials are enabled.
     */
    public void setEnabled(boolean enabled) {
        this.enabled = enabled;
    }

    /**
     * Indicates that the specified potential should not contribute to potential
     * calculations. If potential is not in this group, no action is taken.
     */
    public void setEnabled(IPotentialMolecular potential, boolean enabled) {
        for (PotentialLinker link : potentialList) {
            if (link.potential == potential) {
                link.enabled = enabled;
                return;
            }
        }
    }

    /**
     * Indicates that the specified potential should not contribute to potential
     * calculations. If potential is not in this group, no action is taken.
     */
    public void setEnabled(IPotentialAtomic potential, boolean enabled) {
        for (PotentialLinker link : potentialList) {
            if (link.potential instanceof PotentialGroup && ((PotentialGroup) link.potential).contains(potential)) {
                ((PotentialGroup) link.potential).setEnabled(potential, enabled);
                return;
            }
        }
    }

    /**
     * Returns true if the potential is in this group and has not been disabled
     * via a previous call to setEnabled; returns false otherwise.
     */
    public boolean isEnabled(IPotentialAtomic potential) {
        for (PotentialLinker link : potentialList) {
            if (link.potential instanceof PotentialGroup && ((PotentialGroup) link.potential).contains(potential)) {
                return ((PotentialGroup) link.potential).isEnabled(potential);
            }
        }
        throw new RuntimeException("I don't know about " + potential);
    }

    /**
     * Returns true if the potential is in this group and has not been disabled
     * via a previous call to setEnabled; returns false otherwise.
     */
    public boolean isEnabled(IPotentialMolecular potential) {
        for (PotentialLinker link : potentialList) {
            if (link.potential == potential) {
                return link.enabled;
            }
        }
        throw new RuntimeException("I don't know about " + potential);
    }

    /**
     * Returns an array containing all molecular Potentials.
     */
    public IPotentialMolecular[] getPotentials() {
        IPotentialMolecular[] potentials = new IPotentialMolecular[potentialList.size()];
        int i = 0;
        for (PotentialLinker link : potentialList) {
            potentials[i++] = link.potential;
        }
        return potentials;
    }

    public static class MoleculeIterator0 implements MoleculesetIteratorPDT {
        protected final MoleculeArrayList list = new MoleculeArrayList(0);
        protected boolean finished;

        public MoleculeIterator0() {
        }

        public void setBox(Box box) {
        }

        public void setTarget(IMolecule target) {
        }

        public void setDirection(IteratorDirective.Direction direction) {
        }

        public int nBody() {
            return 0;
        }

        public IMoleculeList next() {
            if (finished) {
                return null;
            }
            finished = true;
            return list;
        }

        public void reset() {
            finished = false;
        }

        public int size() {
            return 1;
        }

        public void unset() {
            finished = true;
        }
    }

    //Added
    public static class AtomIterator0 implements AtomsetIteratorBasisDependent {
        protected final AtomArrayList list = new AtomArrayList(0);
        protected boolean finished;

        public AtomIterator0() {
        }

        public void setTarget(IAtom target) {
        }

        public int nBody() {
            return 0;
        }

        public IAtomList next() {
            if (finished) {
                return null;
            }
            finished = true;
            return list;
        }

        public void reset() {
            finished = false;
        }

        public int size() {
            return 1;
        }

        public void unset() {
            finished = true;
        }

        public void setBasis(IMoleculeList atoms) {
        }

        public int basisSize() {
            return 0;
        }

        public boolean haveTarget(IAtom target) {
            return false;
        }
    }


    public static class PotentialLinker {
        public final IPotentialMolecular potential;
        public final MoleculesetIteratorPDT iterator;
        public final ISpecies[] types;
        public boolean enabled = true;

        //Constructors
        public PotentialLinker(IPotentialMolecular a, MoleculesetIteratorPDT i, ISpecies[] t) {
            potential = a;
            iterator = i;
            if (t != null) {
                types = t.clone();
            } else {
                types = null;
            }
        }
    }//end of PotentialLinker
}//end of PotentialMaster
    
