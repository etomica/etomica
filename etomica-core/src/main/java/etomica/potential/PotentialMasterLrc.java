/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.AtomArrayList;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.molecule.MoleculeArrayList;

import java.util.ArrayList;

/**
 * Collects potentials used for long-range correction.
 * One instance of this class is added to the PotentialMaster of
 * a Simulation when the first Potential0Lrc class is instantiated.
 * All subsequently created Potential0Lrc classes are added to this instance.
 *
 * @see Potential0Lrc
 * @author David Kofke
 */
 
public class PotentialMasterLrc {

    protected PotentialMasterLrc() {
        atomicPotentials = new ArrayList<IPotentialAtomic>();
        molecularPotentials = new ArrayList<IPotentialMolecular>();
        list0 = new MoleculeArrayList(0);
        listLeaf0 = new AtomArrayList(0);
    }
    
    /**
     * Performs given PotentialCalculation on all LRC potentials added to this group.
     * Checks that group is enabled, box is not null, that it has lrcEnabled,
     * and that the given IteratorDirective has includeLrc set to true; if all
     * are so, calculation is performed.
     */
    public void calculate(Box box, IteratorDirective id, PotentialCalculation pc) {
        if(!enabled || !id.includeLrc) return;
        IAtom targetAtom = id.getTargetAtom();
        IMolecule targetMolecule = id.getTargetMolecule();

        if (pc instanceof PotentialCalculationMolecular && (targetAtom != null)) {
            int numMolecularPotentials = molecularPotentials.size();
            for (int i=0; i<numMolecularPotentials; i++) {
                final IPotentialMolecular potential = molecularPotentials.get(i);
                potential.setBox(box);
                ((IPotential0MoleculeLrc)potential).setTargetMolecule(targetMolecule);
                ((PotentialCalculationMolecular)pc).doCalculation(list0, potential);
            }
        }
        int numAtomicPotentials = atomicPotentials.size();
        for (int i=0; i<numAtomicPotentials; i++) {
            final IPotentialAtomic potential = atomicPotentials.get(i);
            potential.setBox(box);
            if (targetMolecule != null) {
                ((IPotential0Lrc)potential).setTargetMolecule(targetMolecule);
            }
            else {
                ((IPotential0Lrc)potential).setTargetAtom(targetAtom);
            }
            pc.doCalculation(listLeaf0, potential);
        }
    }
    
    public boolean isEnabled() {
        return enabled;
    }

    public void setEnabled(boolean enabled) {
        this.enabled = enabled;
    }
    
    public void addPotential(IPotentialMolecular potential) {
        molecularPotentials.add(potential);
    }

    public void addPotential(IPotentialAtomic potential) {
        atomicPotentials.add(potential);
    }

    protected final ArrayList<IPotentialMolecular> molecularPotentials;
    protected final ArrayList<IPotentialAtomic> atomicPotentials;
    protected boolean enabled = true;
    protected final MoleculeArrayList list0;
    protected final AtomArrayList listLeaf0;
    
    private static final long serialVersionUID = 1L;
}
