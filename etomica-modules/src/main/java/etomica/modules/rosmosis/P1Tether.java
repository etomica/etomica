/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.rosmosis;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.box.Box;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.potential.Potential1;
import etomica.potential.PotentialSoft;
import etomica.space.Space;
import etomica.space.Tensor;


/**
 * Potential for a harmonic tether.
 * 
 * U = 0.5 eps r^2
 * 
 * @author Andrew Schultz
 */
public class P1Tether extends Potential1 implements AgentSource<IVector>, PotentialSoft {

    public P1Tether(Box box, ISpecies species, Space _space) {
        super(_space);
        this.species = species;
        agentManager = new AtomLeafAgentManager<IVector>(this, box, IVector.class);
        work = _space.makeVector();
        gradient = new IVector[]{work};
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

    public IVector[] gradient(IAtomList atoms) {
        IAtom atom = atoms.getAtom(0);
        work.E(atom.getPosition());
        work.ME(agentManager.getAgent(atom));
        work.TE(epsilon);
        return gradient;
    }

    public IVector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        return gradient(atoms);
    }

    public double virial(IAtomList atoms) {
        return 0;
    }
    public IVector makeAgent(IAtom a, Box agentBox) {
        if (a.getType().getSpecies() == species) {
            IVector vec = space.makeVector();
            vec.E(a.getPosition());
            return vec;
        }
        return null;
    }

    public void releaseAgent(IVector agent, IAtom atom, Box agentBox) {
        /* do nothing */
    }

    protected final AtomLeafAgentManager<IVector> agentManager;
    protected final ISpecies species;
    protected final IVector work;
    protected final IVector[] gradient;
    protected double epsilon;
}
