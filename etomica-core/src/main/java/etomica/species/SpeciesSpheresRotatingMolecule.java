/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.species;

import etomica.api.IMolecule;
import etomica.atom.*;
import etomica.chem.elements.ElementSimple;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Species in which molecules are made of a single atom.  The molecule itself
 * holds the orientation.
 *
 * @author Andrew Schultz
 */
public class SpeciesSpheresRotatingMolecule extends SpeciesSpheresMono implements ISpeciesOriented {

    private static final long serialVersionUID = 1L;
    protected Vector moment;

    public SpeciesSpheresRotatingMolecule(Simulation sim, Space _space) {
        this(sim, _space, makeNominalMoment(_space));
    }

    public SpeciesSpheresRotatingMolecule(Simulation sim, Space _space, Vector moment) {
        this(_space, new AtomType(new ElementSimple(sim)), moment);
    }

    public SpeciesSpheresRotatingMolecule(Space _space, AtomType atomType, Vector moment) {
        super(_space, atomType);
        this.moment = _space.makeVector();
        this.moment.E(moment);
    }

    protected static final Vector makeNominalMoment(Space space) {
        Vector m = space.makeVector();
        m.E(1);
        return m;
    }

    /**
     * Constructs a new group.
     */
     public IMolecule makeMolecule() {
         MoleculeOriented group = isDynamic ? new MoleculeOrientedDynamic(space, this, 1)
                                            : new MoleculeOriented(space, this, 1);
         group.addChildAtom(makeLeafAtom());
         return group;
     }

    protected IAtom makeLeafAtom() {
        return new Atom(space, leafAtomType);
    }

    public double getMass() {
        return leafAtomType.getMass();
    }

    public Vector getMomentOfInertia() {
        return moment;
    }

    /**
     * Sets the species' moment of inertia to the given moment.
     */
    public void setMomentOfInertia(Vector newMoment) {
        moment.E(newMoment);
    }
}
