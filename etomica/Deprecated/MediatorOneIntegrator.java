package simulate;

/**
* Element coordinator appropriate to use of a single integator and a single phase.
* OK for multiple phases if integrator is appropriate, and if only one phase is displayed.
* Registers all displays and meters as interval listeners of the integrator.
* Registers all phases with the integrator.
* Adds to phase any meter that has phase == null.
* Sets phase for all displays to be the first instantiated phase.
* Adds integrator to the controller.
* This is the default ElementCoordinator.
*/
public class MediatorOneIntegrator extends MediatorBasic {
    
    public MediatorOneIntegrator() {this(Simulation.instance);}
    public MediatorOneIntegrator(Simulation sim) {super(sim);}
    
    public void go() {
        if(completed) return;
        super.go();
        Integrator integrator = parentSimulation().integrator(0);
        //Connect phases and integrator
                //if(phase.integrator() == null && Integrator.this.wantsPhase()) {  //phase doesn't have an integrator and integrator wants a phase
        for(java.util.Iterator i=parentSimulation().phaseList.iterator(); i.hasNext(); ) {
            Phase phase = (Phase)i.next();
            phase.setIntegrator(integrator);
        }
                    
        //Register displays as listeners for integrator interval events and set their phase
        Phase phase = (Phase)parentSimulation().phase(0);
        for(java.util.Iterator i=parentSimulation().displayList.iterator(); i.hasNext(); ) {
            Display display = (Display)i.next();
            integrator.addIntervalListener(display);
            display.setPhase(phase);
        }
                    
        //Set phase of meters, which also registers them as listeners for integrator interval events 
        for(java.util.Iterator i=parentSimulation().meterList.iterator(); i.hasNext(); ) {
            MeterAbstract meter = (MeterAbstract)i.next();
            if(meter.getPhase() == null) meter.setPhase(phase);
        }
                    
        //Add integrator to the controller (assume there is only one since there is only one integrator)
        ((Controller)parentSimulation().controllerList.getFirst()).add(integrator);
                    
    }//end of go method
}//end of CoordinatorOneIntegrator
