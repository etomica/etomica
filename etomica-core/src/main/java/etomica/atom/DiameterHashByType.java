/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.simulation.Simulation;

/**
 * This class hashes atomic diameters based on the atom type.
 * 
 * @author Andrew Schultz
 */
public class DiameterHashByType implements DiameterHash, AtomTypeAgentManager.AgentSource {

    protected final AtomTypeAgentManager agentManager;
    
    public DiameterHashByType(Simulation sim) {
        agentManager = new AtomTypeAgentManager(this, sim);
    }
    
    public double getDiameter(IAtom atom) {
        return getDiameter(atom.getType());
    }

    public double getDiameter(AtomType atomType) {
        return (Double)agentManager.getAgent(atomType);
    }

    public void setDiameter(AtomType type, double newDiameter) {
        agentManager.setAgent(type, newDiameter);
    }

    public Class getSpeciesAgentClass() {
        return Double.class;
    }

    public Object makeAgent(AtomType type) {
        return new Double(-1);
    }

    public void releaseAgent(Object agent, AtomType type) {
    }
}
