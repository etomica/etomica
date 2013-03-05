package etomica.integrator.mcmove;

/**
 * MC move event that indicates the move's trial failed, meaning that
 * MCMove.doTrial() returned false.
 * 
 * @author Andrew Schultz
 */
public class MCMoveTrialFailedEvent extends MCMoveEvent {

    public MCMoveTrialFailedEvent(MCMoveManager moveManager) {
        super();
        this.moveManager = moveManager;
    }
    
    public MCMove getMCMove() {
        return moveManager.getSelectedMove();
    }
    
    private final MCMoveManager moveManager;
}
