package etomica.tmmc;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMove;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;

/**
 * Integrator that implements Transition-Matrix Monte Carlo method.
 *
 * @author David Kofke
 */
 
 /* History of changes
  * 07/10/02 (DAK) started development
  * 09/20/02 (DAK/JS) correction for when doTrial returns false
  */
  
public class IntegratorTMMC extends IntegratorMC {
    
    public IntegratorTMMC(Simulation sim) {
        this(sim.potentialMaster, sim.getDefaults().temperature);
    }
    
    public IntegratorTMMC(PotentialMaster potentialMaster, double temperature) {
        super(potentialMaster, temperature);
        setWeightUpdateInterval(1000000); //10^6
    }
    
    public void setMacrostateManager(MacrostateManager m) {
        macrostateManager = m;
        nStates = macrostateManager.numberOfStates(phase);
        C = new double[nStates][3];
        H = new double[nStates];
        weight = new double[nStates];
    }
    /**
     * Adds a trial that results in a change in the macrostate of the system.
     * Identified as such by adding via this method, rather than simple add
     * method.  Move is treated as any other move, except when invoked additional
     * steps are taken to update transition matrix and weights are used in
     * deciding acceptance.
     */
    /*public void addMacrostateMove(MCMove move) {
        //other stuff
        add(move);
    }*/
    
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
        MCMove move = moveManager.selectMove();
        if(move == null) return;
        
        //perform the trial
        int iStateOld = macrostateManager.stateIndex(phase);//new to tmmc
        if(!move.doTrial()) {
            C[iStateOld][1] += 1;
            if(--doStepCount == 0) updateWeights();
            return;
        }
        
        //notify any listeners that move has been attempted
        if(eventManager != null) { //consider using a final boolean flag that is set in constructor
            event.mcMove = move;
            event.isTrialNotify = true;
            eventManager.fireEvent(event);
        }
        
        //decide acceptance
        int iStateNew = macrostateManager.stateIndex(phase);//new to tmmc
        int iDelta = iStateNew - iStateOld + 1;// 0, 1, 2  new to tmmc
        double weightDifference = weight[iStateNew] - weight[iStateOld]; //new to tmmc
        double lnChi = move.getA() + move.getB();
        double r = (lnChi < 0.0) ? Math.exp(lnChi) : 1.0; //new to tmmc
        C[iStateOld][iDelta] += r;  //new to tmmc
        C[iStateOld][1] += (1.0 - r); //new to tmmc
        lnChi += weightDifference;  //new to tmmc
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
        
        move.getTracker().updateCounts(event.wasAccepted,r);
        
        if(--doStepCount == 0) updateWeights();
    }//end of doStep
    
    private void updateWeights() {
        for(int i=0; i<nStates; i++) {
            H[i] = C[i][0] + C[i][1] + C[i][2];
        }
        weight[0] = 0.0;
        for(int i=1; i<nStates; i++) {
            //w_i = w_(i-1) - log( C[(N-1)->N] / C[N->(N-1)]
            weight[i] = weight[i-1] - Math.log((C[i-1][2]/H[i-1])/(C[i][0]/H[i]));
        }
        doStepCount = weightUpdateInterval;
    }//end of updateWeights 
    
    /**
     * Sets the number of doStep calls between updating of weights.
     */
    public void setWeightUpdateInterval(int i) {
        weightUpdateInterval = i;
        if(weightUpdateInterval < 1) weightUpdateInterval = 1;
        doStepCount = weightUpdateInterval;
    }
    /**
     * Accessor method for number of doStep calls between updating of weights.
     */
    public int getWeightUpdateInterval() {return weightUpdateInterval;}
    
    private double[][] C;
    private double[] H;
    protected double[] weight;
    private int weightUpdateInterval;
    private int doStepCount;
    private MacrostateManager macrostateManager;
    protected int nStates;
}