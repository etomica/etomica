package simulate;
import java.util.*;

public class IntegrationIntervalEvent extends EventObject{
    
    public PhaseSpace phaseSpace;
    public double time;
    public Integrator integrator;
    
    public IntegrationIntervalEvent(Integrator integrator, PhaseSpace phaseSpace) {
        this(integrator,phaseSpace,0.0);
    }
    
    public IntegrationIntervalEvent(Integrator integrator, PhaseSpace phaseSpace, double time) {
        super(integrator);
        this.integrator = integrator;
        this.phaseSpace = phaseSpace;
        this.time = time;
    }
    
}
