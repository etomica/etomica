package etomica.action;

import etomica.*;

public final class SimulationRestart extends SimulationAction {
    
    public SimulationRestart() {
        setLabel("Restart");
    }
    public SimulationRestart(Simulation sim) {
        this();
        setSimulation(sim);
    }
    
    public static void doAction(Simulation sim) {
        for(java.util.Iterator iter=sim.phaseList.iterator(); iter.hasNext(); ) {
            Phase phase = (Phase)iter.next();
            phase.getConfiguration().initializeCoordinates(phase);
        }
        for(java.util.Iterator iter=sim.integratorList.iterator(); iter.hasNext(); ) {
            Integrator integrator = (Integrator)iter.next();
            integrator.halt();//request integrator to stop
            integrator.join();//wait till it does
            integrator.reset();
        }
        for(java.util.Iterator iter=sim.controllerList.iterator(); iter.hasNext(); ) {
            Controller controller = (Controller)iter.next();
            controller.reset();
        }
        
        for(java.util.Iterator iter=sim.meterList.iterator(); iter.hasNext(); ) {
            MeterAbstract meter = (MeterAbstract)iter.next();
            meter.reset();
        }
        for(java.util.Iterator iter=sim.displayList.iterator(); iter.hasNext(); ) {
            Display display = (Display)iter.next();
            display.repaint();
        }
    }
    
    public void actionPerformed(Simulation sim) {
        SimulationRestart.doAction(sim);
    }
}
        