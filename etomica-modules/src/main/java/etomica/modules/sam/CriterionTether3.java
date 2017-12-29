/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.sam;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeAgentManager;
import etomica.molecule.MoleculeAgentManager.MoleculeAgentSource;
import etomica.nbr.NeighborCriterion;
import etomica.simulation.Simulation;
import etomica.species.ISpecies;

/**
 * Returns first leaf atom of each polymer molecule and the atom its bonded to.
 */
public class CriterionTether3 implements NeighborCriterion, MoleculeAgentSource {

    protected final Simulation sim;
    protected final ISpecies polymerSpecies;
    protected final AtomType surfaceType;
    protected Box box;
    protected IMoleculeList polymerList;
    protected MoleculeAgentManager bondManager;
    protected int cursor;
    protected int surfaceCursor;
    protected IMolecule targetMolecule;

    public CriterionTether3(Simulation sim, ISpecies polymerSpecies, AtomType surfaceType) {
        this.sim = sim;
        this.polymerSpecies = polymerSpecies;
        this.surfaceType = surfaceType;
    }

    public void setBox(Box newBox) {
        if (box != newBox) {
            if (box != null) {
                bondManager.dispose();
            }
            bondManager = new MoleculeAgentManager(sim, newBox, this);
        }
        box = newBox;
        polymerList = box.getMoleculeList(polymerSpecies);
    }

    public void setBondedSurfaceAtoms(IMolecule polymerMolecule, IAtomList surfaceAtoms) {
        bondManager.setAgent(polymerMolecule, surfaceAtoms);
    }

    public boolean accept(IAtomList pair) {
        IAtom atom1 = pair.getAtom(0);
        IAtom atom2 = pair.getAtom(1);
        if (atom1.getIndex() != 0 || atom1.getParentGroup().getType() != polymerSpecies) {
            IAtom foo = atom2;
            atom2 = atom1;
            atom1 = foo;
            if (atom1.getIndex() != 0 || atom1.getParentGroup().getType() != polymerSpecies) {
                return false;
            }
        }
        if (atom2.getType() != surfaceType) {
            return false;
        }
        IAtomList bondedSurfaceAtoms = ((IAtomList)bondManager.getAgent(atom1.getParentGroup()));
        if (bondedSurfaceAtoms == null) {
            return false;
        }
        for (int i=0; i<bondedSurfaceAtoms.getAtomCount(); i++) {
            if (bondedSurfaceAtoms.getAtom(i) == atom2) {
                return true;
            }
        }
        return false;
    }

    public boolean needUpdate(IAtom atom) {
        return false;
    }

    public void reset(IAtom atom) {
    }

    public boolean unsafe() {
        return false;
    }

    public Object makeAgent(IMolecule a) {
        return null;
    }

    public void releaseAgent(Object agent, IMolecule atom) {
    }
}
