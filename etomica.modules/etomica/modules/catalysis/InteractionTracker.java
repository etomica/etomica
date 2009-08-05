package etomica.modules.catalysis;

import etomica.api.IAtom;
import etomica.api.IBox;
import etomica.api.ISpecies;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.integrator.IntegratorHard.Agent;
import etomica.integrator.IntegratorHard.CollisionListener;

public class InteractionTracker implements CollisionListener, AgentSource {

    public InteractionTracker(IBox box, ISpecies speciesSurface) {
        agentManager = new AtomLeafAgentManager(this, box);
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
            ((CatalysisAgent)agentManager.getAgent(gasAtom)).nSurfaceBonds++;
        }
        else if (de > 0) {
            // escape
            ((CatalysisAgent)agentManager.getAgent(gasAtom)).nSurfaceBonds--;
        }
    }
    
    public AtomLeafAgentManager getAgentManager() {
        return agentManager;
    }

    public Class getAgentClass() {
        return CatalysisAgent.class;
    }

    public Object makeAgent(IAtom a) {
        if (a.getType().getSpecies() == speciesSurface) return null;
        return new CatalysisAgent();
    }

    public void releaseAgent(Object agent, IAtom atom) {}

    protected final AtomLeafAgentManager agentManager;
    protected final ISpecies speciesSurface;

    public static class CatalysisAgent {
        public int nSurfaceBonds;
        public boolean isRadical;
        public IAtom bondedAtom1;
        public IAtom bondedAtom2;
    }
}
