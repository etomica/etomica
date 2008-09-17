package etomica.integrator.mcmove;

import etomica.api.IEvent;

public abstract class MCMoveEvent implements IEvent {
    
    public abstract MCMove getMCMove();
    
}