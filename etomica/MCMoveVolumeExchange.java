package etomica;

import java.util.Random;

public final class MCMoveVolumeExchange extends MCMove {
    
    private final Random rand = new Random();
    private Phase firstPhase;
    private Phase secondPhase;
    private PhaseAction.Inflate inflate1;
    private PhaseAction.Inflate inflate2;
    private final double ROOT;

    public MCMoveVolumeExchange(IntegratorMC parent) {
        super();
        parentIntegrator = parent;
        ROOT = 1.0/(double)parentIntegrator.parentSimulation().space().D();
        setStepSizeMax(Double.MAX_VALUE);
        setStepSizeMin(Double.MIN_VALUE);
        setStepSize(0.3);
    }
    
    /**
     * Overrides superclass method so that it performs no action.
     * Must set using method that takes an array of phases.
     */
    public void setPhase(Phase p) {}

    public void setPhase(Phase[] p) {
        super.setPhase(p);
        firstPhase = p[0];
        secondPhase = p[1];
        if(firstPhase == null && secondPhase == null) return;
        inflate1 = new PhaseAction.Inflate(firstPhase);
        inflate2 = new PhaseAction.Inflate(secondPhase);
    }
    
    //under revision--- does not work for multiatomics, since intramolecular energy is not considered
    public void thisTrial() {
        double hOld = firstPhase.energy.potential() + secondPhase.energy.potential();
        double v1Old = firstPhase.volume();
        double v2Old = secondPhase.volume();
        double vRatio = v1Old/v2Old * Math.exp(stepSize*(rand.nextDouble() - 0.5));
        double v2New = (v1Old + v2Old)/(1 + vRatio);
        double v1New = (v1Old + v2Old - v2New);
        double v1Scale = v1New/v1Old;
        double v2Scale = v2New/v2Old;
//        System.out.println(v1Scale + "  " + v2Scale + " " +stepSize + " "+nTrials+" "+nAccept+" "+adjustInterval+" "+frequency);
        inflate1.setScale(Math.pow(v1Scale,ROOT));
        inflate2.setScale(Math.pow(v2Scale,ROOT));
        inflate1.attempt();
        inflate2.attempt();
        double hNew = firstPhase.energy.potential() + secondPhase.energy.potential();
        if(hNew >= Double.MAX_VALUE ||
             Math.exp(-(hNew-hOld)/parentIntegrator.temperature+
                       (firstPhase.moleculeCount+1)*Math.log(v1Scale) +
                       (secondPhase.moleculeCount+1)*Math.log(v2Scale))
                < rand.nextDouble()) 
            {  //reject
              inflate1.undo();
              inflate2.undo();
        }
        else nAccept++;
    }//end of thisTrial
}//end of MCMoveVolumeExchange