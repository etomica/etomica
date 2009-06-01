package etomica.space;

import etomica.api.IBoundary;
import etomica.api.IBoundaryEvent;

public class BoundaryEvent implements IBoundaryEvent {

    protected IBoundary boundary = null;
    
    public BoundaryEvent(IBoundary _boundary) {
        boundary = _boundary;
    }
    
    public IBoundary getBoundary() {
        return boundary;
    }
}
