/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.catalysis;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.modules.catalysis.InteractionTracker.CatalysisAgent;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.units.Kelvin;
import etomica.util.random.IRandom;

public class ReactionManagerCO implements IntegratorListener {

    public ReactionManagerCO(Catalysis sim) {
        this.sim = sim;
    }
    
    public void integratorInitialized(IntegratorEvent e) {
    }

    public void integratorStepFinished(IntegratorEvent e) {
        Box box = sim.getBox(0);
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
                if (random.nextDouble() < sim.integrator.getTimeStep()*Math.exp(-uReactCORev/temperature)) {
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

    public void integratorStepStarted(IntegratorEvent e) {
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

    public double getuReactCORev() {
        return uReactCORev;
    }

    public void setuReactCORev(double newUReactCORev) {
        uReactCORev = newUReactCORev;
    }

    protected final Catalysis sim;
    protected int nReactCO = 3;      // number of surface sites C needs for radicalization
    protected double uReactCO = Kelvin.UNIT.toSim(80); // barrier for forward reaction (radicalization)
    protected double uReactCORev = Kelvin.UNIT.toSim(40); // barrier for reverse reaction
}
