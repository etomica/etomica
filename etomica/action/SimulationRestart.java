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
        
        for(Iterator iter=sim.phaseList().iterator(); iter.hasNext(); ) {
            Phase phase = (Phase)iter.next();
            phase.getConfiguration().initializeCoordinates(phase.speciesMaster().childAtomArray());
        }
        for(Iterator iter=sim.integratorList().iterator(); iter.hasNext(); ) {
            Integrator integrator = (Integrator)iter.next();
            integrator.halt();//request integrator to stop
            integrator.join();//wait till it does
            integrator.reset();
        }
        for(Iterator iter=sim.controllerList().iterator(); iter.hasNext(); ) {
            Controller controller = (Controller)iter.next();
            controller.reset();
        }
        
        for(Iterator iter=sim.meterList().iterator(); iter.hasNext(); ) {
            MeterAbstract meter = (MeterAbstract)iter.next();
            //take care that histxxx is reset but restored to prior OnMeterReset condition
            boolean resetHistory = meter.isResetHistoryOnMeterReset();
            boolean resetHistogram = meter.isResetHistogramOnMeterReset();
            meter.setResetHistoryOnMeterReset(true);
            meter.setResetHistogramOnMeterReset(true);
            meter.reset();
            meter.setResetHistoryOnMeterReset(resetHistory);
            meter.setResetHistogramOnMeterReset(resetHistogram);
        }
/*        for(Iterator iter=sim.displayList().iterator(); iter.hasNext(); ) {
            Display display = (Display)iter.next();
            display.doUpdate();
            display.graphic().repaint();
        }*/
    }
    
    public void actionPerformed(Simulation sim) {
        SimulationRestart.doAction(sim);
    }
}
        