package etomica;

import java.util.Random;

public abstract class MCMove implements Cloneable, java.io.Serializable {
    
    int frequency, nominalFrequency;
    double acceptanceRatio, acceptanceTarget, stepSize, stepSizeMax, stepSizeMin;
    int nTrials, nAccept, nTrialsSum, nAcceptSum, adjustInterval;
    boolean perParticleFrequency;
    IntegratorMC parentIntegrator;
    protected MCMove nextMove;
    protected boolean tunable = true;
    protected Phase phase;
    private String name;
    
    public MCMove() {
        nTrials = 0;
        nAccept = 0;
        setAcceptanceTarget(0.5);
        setFrequency(1);
        setPerParticleFrequency(false);
        setAdjustInterval(100);
    }
    
    public void doTrial() {
        nTrials++;
        thisTrial();
        if(nTrials > adjustInterval*frequency) {adjustStepSize();}
    }
    
    public void setParentIntegrator(IntegratorMC parent) {parentIntegrator = parent;}
    public IntegratorMC parentIntegrator() {return parentIntegrator;}
    
    public abstract void thisTrial();
    
    public void setPhase(Phase[] p) {
        setPhase(p[0]);
    }
    
    public void setPhase(Phase p) {
        phase = p;
    }
    
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
    
    
    public void setFrequency(int f) {
        nominalFrequency = f;
        frequency = f;
    }
    public final int getFrequency() {return frequency;}
    public void resetFrequency() {
        frequency = (perParticleFrequency && phase!=null) ? nominalFrequency*phase.moleculeCount : nominalFrequency;
    }
    public final void setPerParticleFrequency(boolean b) {perParticleFrequency = b;}
    public final boolean isPerParticleFrequency() {return perParticleFrequency;}
    
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
    
    public Object clone() {
        Object o = null;
        try {
            o = super.clone();
        } catch(CloneNotSupportedException e) {}
        return o;
    }
    
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