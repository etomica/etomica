package simulate;
import java.util.*;

public class IntegrationIntervalEvent extends EventObject{
    
    public Phase phase;
    public double time;
    public Integrator integrator;
    
    public IntegrationIntervalEvent(Integrator integrator, Phase phase) {
        this(integrator,phase,0.0);
    }
    
    public IntegrationIntervalEvent(Integrator integrator, Phase phase, double time) {
        super(integrator);
        this.integrator = integrator;
        this.phase = phase;
        this.time = time;
    }
    
}
