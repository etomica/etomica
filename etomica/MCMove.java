package etomica;

public abstract class MCMove implements java.io.Serializable {
    
    private int frequency, fullFrequency;
    double acceptanceRatio, acceptanceTarget;
    protected double stepSize, stepSizeMax, stepSizeMin;
    private int nTrials, nAccept, nTrialsSum, nAcceptSum, adjustInterval;
    boolean perParticleFrequency;
    protected MCMove nextMove;
    protected boolean tunable = true;
    protected Phase phase;//should get rid of this
    private Phase[] phases;
    private String name;
    protected final IntegratorMC parentIntegrator;
    protected final PotentialMaster potential;
    protected final PotentialCalculationEnergySum energy;
    
    public MCMove(IntegratorMC parent) {
        parentIntegrator = parent;
        potential = parentIntegrator.potential;
        energy = parent.parentSimulation().energySum(this);
        parent.add(this);
        nTrials = 0;
        nAccept = 0;
        setAcceptanceTarget(0.5);
        setFrequency(1);
        setPerParticleFrequency(false);
        setAdjustInterval(100);
    }
    
    /**
     * Updates statistics regarding the rate of acceptance of trials
     * for this move.
     */
    public void updateCounts(boolean moveWasAccepted) {
        nTrials++;
        if(moveWasAccepted) nAccept++;
        if(nTrials > adjustInterval*frequency) {adjustStepSize();}
    }
    
    public IntegratorMC parentIntegrator() {return parentIntegrator;}
    
    /**
     * Method to perform trial move.  Returns false if the trial could
     * not be attempted, for example if there were no molecules in the
     * phase and the trial is designed to displace an atom; returns true otherwise.
     */
    public abstract boolean doTrial();
    
    /**
     * Returns log of the ratio of the trial probabilities, ln(Tji/Tij) for the
     * states encountered before (i) and after (j) the most recent call to doTrial(). 
     * Tij is the probability that this move would generate state j from state i, and
     * Tji is the probability that a subsequent call to doTrial would return to state i
     * from state j.
     */
    public abstract double lnTrialRatio();
    
    /**
     * Returns the log of the limiting-distribution probabilities of states, ln(Pj/Pi), 
     * for the states encountered before (i) and after (j) the most recent call to 
     * doTrial.
     */
    public abstract double lnProbabilityRatio();
    
    /**
     * Method called by IntegratorMC in the event that the most recent trial is accepted.
     */
    public abstract void acceptNotify();
    
    /**
     * Method called by IntegratorMC in the event that the most recent trial move is
     * rejected.  This method should cause the system to be restored to the condition
     * before the most recent call to doTrial.
     */
    public abstract void rejectNotify();
    
    /**
     * Returns an iterator that yields the atoms that were affected by
     * the trial move the last time thisTrial was invoked 
     * (regardless of whether the move was accepted).  This information usually
     * is not needed, but it is available in cases where required by objects that
     * perform some activity while the move is in progress, or need to update 
     * after the move is completed.  Such objects can receive notification of the
     * move's completion  by registering with the IntegratorMC as MCMoveEventListeners.
     */
    public abstract AtomIterator affectedAtoms();
    
    public void setPhase(Phase[] p) {
        phases = p;
        setPhase(p[0]);
    }
    
    public void setPhase(Phase p) {
        phase = p;
        if(phases[0] != p) {
            phases = new Phase[] {p};
        }
    }
    public Phase[] getPhases() {return phases;}
    
    public void adjustStepSize() {
        if(nTrials == 0) {return;}
        nTrialsSum += nTrials;
        nAcceptSum += nAccept;
        if(nTrialsSum != 0) acceptanceRatio = (double)nAcceptSum/(double)nTrialsSum;
        if(!tunable) return;
        if(nAccept > (int)(acceptanceTarget*nTrials)) {stepSize *= 1.05;}
        else{stepSize *= 0.95;}
        stepSize = Math.min(stepSize,stepSizeMax);
        stepSize = Math.max(stepSize,stepSizeMin);
        nTrials = 0;
        nAccept = 0;
    }
    
    /**
     * Set an unnormalized frequency for performing this move, relative to the
     * other moves that have been added to the integrator.  Each move is performed
     * (on average) an amount in proportion to this frequency.  Moves having
     * the same frequency are performed with equal likelihood.  Frequency may
     * be adjusted as specified by the perParticleFrequency flag.
     * @see #setPerParticleFrequency
     */
    public void setFrequency(int f) {
        frequency = f;
        if(parentIntegrator != null) parentIntegrator.doReset();
    }
    /**
     * Accessor method for move frequency.
     * @see #setFrequency
     */
    public final int getFrequency() {return frequency;}
    
    /**
     * Indicates if frequency indicates total frequency, or frequency
     * per particle.  If per particle, frequency input to setFrequency
     * method is multiplied by the number of particles affected by the move
     * (not regularly updated for changing particle numbers).
     * @see #setFrequency
     */ 
    public final void setPerParticleFrequency(boolean b) {
        perParticleFrequency = b;
        if(parentIntegrator != null) parentIntegrator.doReset();
    }
    /**
     * Accessor method for perParticleFrequency flag.
     * @see #setPerParticleFrequency
     */
    public final boolean isPerParticleFrequency() {return perParticleFrequency;}
    
    /**
     * Frequency this move is performed, considering the perParticleFrequency flag.
     * Used by IntegratorMC
     */
    int fullFrequency() {return fullFrequency;}
    
    /**
     * Updates the full frequency based on the current value of the frequency,
     * the status of the perParticleFrequency flag, and the current number of
     * molecules in the phases affected by the move.
     * Invoked by IntegratorMC.doReset
     */
    void resetFullFrequency() {
        fullFrequency = frequency;
        if(perParticleFrequency && phases!=null) {
            int mCount = 0;
            for(int i=0; i<phases.length; i++) 
                if(phases[i]!= null) mCount += phases[i].moleculeCount();
            fullFrequency *= mCount;
        }
    }
            
    
    /**
     * Fraction of time trials of this type were accepted since acceptanceTarget was set
     */
    public double acceptanceRatio() {return acceptanceRatio;}
    public final void setAcceptanceTarget(double a) {
        nTrialsSum = 0;
        nAcceptSum = 0;
        acceptanceTarget = a;
    }
    public final double getAcceptanceTarget() {return acceptanceTarget;}
    
    public final void setStepSize(double step) {stepSize = step;}
    public final double getStepSize() {return stepSize;}
    
    public final void setStepSizeMax(double step) {stepSizeMax = step;}
    public final double getStepSizeMax() {return stepSizeMax;}
    public final void setStepSizeMin(double step) {stepSizeMin = step;}
    public final double getStepSizeMin() {return stepSizeMin;}
    
    public final void setNextMove(MCMove move) {nextMove = move;}
    public final MCMove nextMove() {return nextMove;}
    
    public final void setAdjustInterval(int i) {adjustInterval = i;}
    public final int getAdjustInterval() {return adjustInterval;}
    
    /**
     * Sets a flag to indicate whether tuning of the move is to be performed
     * Tuning aims to set the acceptance rate to some target value
     * Some moves (e.g., simple insertion trial) are inherently untunable
     */
    public final void setTunable(boolean b) {tunable = b;}
    public final boolean getTunable() {return tunable;}
    
    /**
     * Accessor method of the name of this object
     * 
     * @return The given name
     */
    public final String getName() {return name;}

    /**
     * Method to set the name of this object
     * 
     * @param name The name string to be associated with this object
     */
    public final void setName(String name) {this.name = name;}

    /**
     * Overrides the Object class toString method to have it return the output of getName
     * 
     * @return The name given to the object
     */
    public String toString() {return getName();}  //override Object method
          
    
}