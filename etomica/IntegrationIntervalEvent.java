package simulate;
import java.util.*;

public class IntegrationIntervalEvent extends EventObject{
    
    Phase phase;
    double time;
    
    public IntegrationIntervalEvent(Integrator integrator, Phase phase) {
        super(integrator);
        this.phase = phase;
        time = 0.0;
    }
    
    public IntegrationIntervalEvent(Integrator integrator, Phase phase, double time) {
        super(integrator);
        this.phase = phase;
        this.time = time;
    }
    
}
