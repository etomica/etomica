/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.rosmosis;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.potential.IPotentialField;
import etomica.potential.Potential1;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;


/**
 * Potential for a harmonic tether.
 * <p>
 * U = 0.5 eps r^2
 *
 * @author Andrew Schultz
 */
public class P1Tether extends Potential1 implements AgentSource<Vector>, IPotentialField {

    public P1Tether(Box box, ISpecies species, Space _space) {
        super(_space);
        this.space = _space;
        this.species = species;
        agentManager = new AtomLeafAgentManager<Vector>(this, box);
        work = _space.makeVector();
    }

    public void setEpsilon(double newEpsilon) {
        epsilon = newEpsilon;
    }

    public double getEpsilon() {
        return epsilon;
    }

    public double u(IAtom atom) {
        return 0.5 * epsilon * atom.getPosition().Mv1Squared(agentManager.getAgent(atom));
    }

    public double energy(IAtomList atoms) {
        return u(atoms.get(0));
    }

    public double udu(IAtom atom, Vector f) {
        work.E(atom.getPosition());
        work.ME(agentManager.getAgent(atom));
        double x2 = work.squared();
        work.TE(epsilon);
        f.ME(work);
        return 0.5 * epsilon * x2;
    }

    public Vector makeAgent(IAtom a, Box agentBox) {
        if (a.getType().getSpecies() == species) {
            Vector vec = space.makeVector();
            vec.E(a.getPosition());
            return vec;
        }
        return null;
    }

    public void releaseAgent(Vector agent, IAtom atom, Box agentBox) {
        /* do nothing */
    }

    protected final Space space;
    protected final AtomLeafAgentManager<Vector> agentManager;
    protected final ISpecies species;
    protected final Vector work;
    protected double epsilon;
}
