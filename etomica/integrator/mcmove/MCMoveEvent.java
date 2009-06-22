package etomica.integrator.mcmove;

import etomica.util.IEvent;

public abstract class MCMoveEvent implements IEvent {
    
    public abstract MCMove getMCMove();
    
}