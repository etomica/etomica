package etomica.action;

import etomica.*;
//Java2 imports
//import java.util.Iterator;

import etomica.utility.Iterator;

public final class SimulationRestart extends SimulationAction {
    
    public SimulationRestart() {
        setLabel("Restart");
    }
    public SimulationRestart(Simulation sim) {
        this();
        setSimulation(sim);
    }
    
    public static void doAction(Simulation sim) {
        for(Iterator iter=sim.phaseList.iterator(); iter.hasNext(); ) {
            Phase phase = (Phase)iter.next();
            phase.getConfiguration().initializeCoordinates(phase);
        }
        for(Iterator iter=sim.integratorList.iterator(); iter.hasNext(); ) {
            Integrator integrator = (Integrator)iter.next();
            integrator.halt();//request integrator to stop
            integrator.join();//wait till it does
            integrator.reset();
        }
        for(Iterator iter=sim.controllerList.iterator(); iter.hasNext(); ) {
            Controller controller = (Controller)iter.next();
            controller.reset();
        }
        
        for(Iterator iter=sim.meterList.iterator(); iter.hasNext(); ) {
            MeterAbstract meter = (MeterAbstract)iter.next();
            meter.reset();
        }
        for(Iterator iter=sim.displayList.iterator(); iter.hasNext(); ) {
            Display display = (Display)iter.next();
            display.repaint();
        }
    }
    
    public void actionPerformed(Simulation sim) {
        SimulationRestart.doAction(sim);
    }
}
        