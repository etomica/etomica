package etomica;

/**
 * Integrator to perform Metropolis Monte Carlo sampling.  Works with a set
 * of MCMove instances that are added to the integrator upon their construction.
 * A step performed by the integrator consists of selecting a MCMove from the
 * set, performing the trial defined by the MCMove, and deciding acceptance
 * of the trial using information from the MCMove.
 *
 * @author David Kofke
 */
public class IntegratorMC extends Integrator implements EtomicaElement {
    
    public String version() {return "IntegratorMC:02.03.19"+Integrator.VERSION;}
    
    private MCMove firstMove, lastMove;
    private int frequencyTotal;
    private int moveCount;
    private SimulationEventManager eventManager;
    private final MCMoveEvent event = new MCMoveEvent(this);
    
    public IntegratorMC() {
        this(Simulation.instance);
    }
    public IntegratorMC(Simulation sim) {
        super(sim);
        setIsothermal(true); //has no practical effect, but sets value of isothermal to be consistent with way integrator is sampling
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("General Monte Carlo simulation");
        return info;
    }
    
    /**
     * Sets moves in given array to be integrator's set of moves, deleting any existing moves.
     */
    public void setMCMoves(MCMove[] moves) {
        firstMove = null;
        moveCount = 0;
        for(int i=0; i<moves.length; i++) {
            add(moves[i]);
        }
    }
    
    /**
     * Constructs and returns array of all moves added to the integrator.  
     */
    public MCMove[] getMCMoves() {
        MCMove[] moves = new MCMove[moveCount];
        int i=0;
        for(MCMove move=firstMove; move!=null; move=move.nextMove()) {moves[i++]=move;}
        return moves;
    }

    /**
     * Adds a basic MCMove to the set of moves performed by the integrator.
     * Called only in constructor of the MCMove class.
     */
    void add(MCMove move) {
        //make sure this is this parent of the move being added
        if(move.parentIntegrator() != this) {
            throw new RuntimeException("Inappropriate move added in IntegratorMC");
        }
        //make sure move wasn't added already
        for(MCMove m=firstMove; m!=null; m=m.nextMove()) {
            if(move == m) {
                throw new RuntimeException("Attempt to add move twice in IntegratorMC; Move is added in its constructor, and should not be added again");
            }
        }
        if(firstMove == null) {firstMove = move;}
        else {lastMove.setNextMove(move);}
        lastMove = move;
        move.setNextMove(null);
        move.setPhase(phase);
        moveCount++;
    }
    
    /**
     * Invokes superclass method and informs all MCMoves about the new phase.
     */
    public boolean addPhase(Phase p) {
        if(!super.addPhase(p)) return false;
        for(MCMove move=firstMove; move!=null; move=move.nextMove()) {move.setPhase(phase);}
        return true;
    }
    
    /**
     * Method to select and perform an elementary Monte Carlo move.  
     * The type of move performed is chosen from all MCMoves that have been added to the
     * integrator.  Each MCMove has associated with it a (unnormalized) frequency, which
     * when weighed against the frequencies given the other MCMoves, determines
     * the likelihood that the move is selected.
     * After completing move, fires an MCMove event if there are any listeners.
     */
    public void doStep() {
        //select the move
        if(firstMove == null) return;
        int i = (int)(Simulation.random.nextDouble()*frequencyTotal);
        MCMove move = firstMove;
        while((i-=move.fullFrequency()) >= 0) {
            move = move.nextMove();
        }
        
        //perform the trial
        if(!move.doTrial()) return;
        
        //notify any listeners that move has been attempted
        if(eventManager != null) { //consider using a final boolean flag that is set in constructor
            event.mcMove = move;
            event.isTrialNotify = true;
            eventManager.fireEvent(event);
        }
        
        //decide acceptance
        double lnChi = move.lnTrialRatio() + move.lnProbabilityRatio();
        if(lnChi <= -Double.MAX_VALUE || 
                (lnChi < 0.0 && Math.exp(lnChi) < Simulation.random.nextDouble())) {//reject
            move.rejectNotify();
            event.wasAccepted = false;
        } else {
            move.acceptNotify();
            event.wasAccepted = true;
        }

        //notify listeners of outcome
        if(eventManager != null) { //consider using a final boolean flag that is set in constructor
            event.isTrialNotify = false;
            eventManager.fireEvent(event);
        }
        
        move.updateCounts(event.wasAccepted);
    }
    
    /**
     * Recomputes all the move frequencies.
     */
    public void doReset() {
        frequencyTotal = 0;
        for(MCMove m=firstMove; m!=null; m=m.nextMove()) {
            m.resetFullFrequency();
            frequencyTotal += m.fullFrequency();
        }
    }
    
    public void addMCMoveListener(MCMoveListener listener) {
        if(eventManager == null) eventManager = new SimulationEventManager();
        eventManager.addListener(listener);
    }
    public void removeMCMoveListener(MCMoveListener listener) {
        if(eventManager == null) return; //should define an exception
        eventManager.removeListener(listener);
    }
    
    public Integrator.Agent makeAgent(Atom a) {
        return null;
    }
}//end of IntegratorMC