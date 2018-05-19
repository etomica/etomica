/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.catalysis;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.IntegratorHard.Agent;
import etomica.integrator.IntegratorHard.CollisionListener;
import etomica.species.ISpecies;

public class InteractionTracker implements CollisionListener, AgentSource<InteractionTracker.CatalysisAgent> {

    public InteractionTracker(Box box, ISpecies speciesSurface) {
        agentManager = new AtomLeafAgentManager<CatalysisAgent>(this, box);
        this.speciesSurface = speciesSurface;
    }
    
    public void collisionAction(Agent colliderAgent) {
        IAtom gasAtom = colliderAgent.collisionPartner;
        if (gasAtom == null) return;
        if (gasAtom.getType().getSpecies() == speciesSurface) {
            // surface atoms don't collide with themselves, so atom must be gas
            gasAtom = colliderAgent.atom;
        }
        else if (colliderAgent.atom.getType().getSpecies() != speciesSurface) {
            // two gas atoms -- we don't care
            return;
        }

        // one of our atoms is a surface atom
        double de = colliderAgent.collisionPotential.energyChange();
        if (de < 0) {
            // capture
            agentManager.getAgent(gasAtom).nSurfaceBonds++;
        }
        else if (de > 0) {
            // escape
            agentManager.getAgent(gasAtom).nSurfaceBonds--;
        }
    }
    
    public void reset() {
        IAtomList list = agentManager.getBox().getLeafList();
        for (int i = 0; i<list.size(); i++) {
            CatalysisAgent agent = agentManager.getAgent(list.get(i));
            if (agent == null) {
                continue;
            }
            agent.nSurfaceBonds = 0;
            agent.isRadical = false;
            // ConfigurationCatalysis is responsible for resetting bonds
        }
    }
    
    public AtomLeafAgentManager<CatalysisAgent> getAgentManager() {
        return agentManager;
    }

    public CatalysisAgent makeAgent(IAtom a, Box agentBox) {
        if (a.getType().getSpecies() == speciesSurface) return null;
        return new CatalysisAgent();
    }

    public void releaseAgent(CatalysisAgent agent, IAtom atom, Box agentBox) {}

    protected final AtomLeafAgentManager<CatalysisAgent> agentManager;
    protected final ISpecies speciesSurface;

    public static class CatalysisAgent {
        public int nSurfaceBonds;
        public boolean isRadical;
        public IAtom bondedAtom1;
        public IAtom bondedAtom2;
    }
}
