package simulate;

import java.awt.Component;
import java.util.Random;

public abstract class MCMove extends Component {
    
    int frequency, nominalFrequency;
    double acceptanceTarget, stepSize, stepSizeMax, stepSizeMin;
    int nTrials, nAccept, nTrialsSum, nAcceptSum, adjustInterval;
    boolean perParticleFrequency;
    Integrator parentIntegrator;
    protected MCMove nextMove;
    
    public MCMove() {
        nTrials = 0;
        nAccept = 0;
        acceptanceTarget = 0.50;
        setFrequency(1);
        setPerParticleFrequency(false);
        setAdjustInterval(100);
    }
    
    public void doTrial(Phase phase) {
        nTrials++;
        thisTrial(phase);
        if(nTrials > adjustInterval*frequency) {adjustStepSize();}
    }
    
    public abstract void thisTrial(Phase phase);
    
    public void adjustStepSize() {
        if(nTrials == 0) {return;}
        if(nAccept > (int)(acceptanceTarget*nTrials)) {stepSize *= 1.05;}
        else{stepSize *= 0.95;}
        stepSize = Math.min(stepSize,stepSizeMax);
        stepSize = Math.max(stepSize,stepSizeMin);
        nTrialsSum += nTrials;
        nAcceptSum += nAccept;
        nTrials = 0;
        nAccept = 0;
    }
    
    public void setFrequency(int f) {
        nominalFrequency = f;
        frequency = f;
    }
    public final int getFrequency() {return frequency;}
    public void resetFrequency(Phase phase) {
        frequency = perParticleFrequency ? nominalFrequency*phase.moleculeCount : nominalFrequency;
    }
    public final void setPerParticleFrequency(boolean b) {perParticleFrequency = b;}
    public final boolean isPerParticleFrequency() {return perParticleFrequency;}
    
    public final void setAcceptanceTarget(double a) {acceptanceTarget = a;}
    public final double getAcceptanceTarget() {return acceptanceTarget;}
    
    public final void setStepSize(double step) {stepSize = step;}
    public final double getStepSize() {return stepSize;}
    
    public final void setStepSizeMax(double step) {stepSizeMax = step;}
    public final double getStepSizeMax() {return stepSizeMax;}
    public final void setStepSizeMin(double step) {stepSizeMin = step;}
    public final double getStepSizeMin() {return stepSizeMin;}
    
    public final void setNextMove(MCMove move) {nextMove = move;}
    public final MCMove getNextMove() {return nextMove;}
    
    public final void setAdjustInterval(int i) {adjustInterval = i;}
    public final int getAdjustInterval() {return adjustInterval;}
    
}