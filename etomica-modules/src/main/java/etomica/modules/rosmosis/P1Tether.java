/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.rosmosis;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.potential.Potential1;
import etomica.potential.PotentialSoft;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.species.ISpecies;


/**
 * Potential for a harmonic tether.
 * 
 * U = 0.5 eps r^2
 * 
 * @author Andrew Schultz
 */
public class P1Tether extends Potential1 implements AgentSource<Vector>, PotentialSoft {

    public P1Tether(Box box, ISpecies species, Space _space) {
        super(_space);
        this.species = species;
        agentManager = new AtomLeafAgentManager<Vector>(this, box, Vector.class);
        work = _space.makeVector();
        gradient = new Vector[]{work};
    }
    
    public void setEpsilon(double newEpsilon) {
        epsilon = newEpsilon;
    }
    
    public double getEpsilon() {
        return epsilon;
    }

    public double energy(IAtomList atoms) {
        IAtom atom = atoms.getAtom(0);
        work.E(atom.getPosition());
        work.ME(agentManager.getAgent(atom));
        return 0.5 * epsilon * work.squared();
    }

    public Vector[] gradient(IAtomList atoms) {
        IAtom atom = atoms.getAtom(0);
        work.E(atom.getPosition());
        work.ME(agentManager.getAgent(atom));
        work.TE(epsilon);
        return gradient;
    }

    public Vector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        return gradient(atoms);
    }

    public double virial(IAtomList atoms) {
        return 0;
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

    protected final AtomLeafAgentManager<Vector> agentManager;
    protected final ISpecies species;
    protected final Vector work;
    protected final Vector[] gradient;
    protected double epsilon;
}
