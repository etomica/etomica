package etomica.modules.catalysis;

import etomica.api.IAtom;
import etomica.api.IBox;
import etomica.api.IIntegratorEvent;
import etomica.api.IIntegratorListener;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.IRandom;
import etomica.atom.AtomLeafAgentManager;
import etomica.modules.catalysis.InteractionTracker.CatalysisAgent;

public class ReactionManagerCO implements IIntegratorListener {

    public ReactionManagerCO(Catalysis sim) {
        this.sim = sim;
    }
    
    public void integratorInitialized(IIntegratorEvent e) {
    }

    public void integratorStepFinished(IIntegratorEvent e) {
        IBox box = sim.getBox(0);
        double temperature = sim.integrator.getTemperature();
        IRandom random = sim.getRandom();
        AtomLeafAgentManager agentManager = sim.interactionTracker.getAgentManager();

        IMoleculeList listC = box.getMoleculeList(sim.speciesC);
        for (int i=0; i<listC.getMoleculeCount(); i++) {
            IMolecule molecule = listC.getMolecule(i);
            IAtom atom = molecule.getChildList().getAtom(0);
            CatalysisAgent agent = (CatalysisAgent)agentManager.getAgent(atom);
            if (agent.bondedAtom1.getType() == atom.getType()) {
                throw new RuntimeException("oops");
            }
            if (agent.isRadical) {
                // consider deradicalization
                if (random.nextDouble() < sim.integrator.getTimeStep()*Math.exp(-uReactCO/temperature)) {
                    agent.isRadical = false;
                }
            }
            else if (agent.nSurfaceBonds >= nReactCO && agent.bondedAtom2 == null) {
                // CO on the surface, consider becoming a radical
                if (random.nextDouble() < sim.integrator.getTimeStep()*Math.exp(-uReactCO/temperature)) {
                    agent.isRadical = true;
                }
            }
        }
   }

    public void integratorStepStarted(IIntegratorEvent e) {
    }

    public int getnReactCO() {
        return nReactCO;
    }

    public void setnReactCO(int newNReactCO) {
        nReactCO = newNReactCO;
    }

    public double getuReactCO() {
        return uReactCO;
    }

    public void setuReactCO(double newUReactCO) {
        uReactCO = newUReactCO;
    }

    protected final Catalysis sim;
    protected int nReactCO = 3;      // number of surface sites C needs for radicalization
    protected double uReactCO = 2.0; // barrier for forward and reverse reactions
}
