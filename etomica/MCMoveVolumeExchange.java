package etomica;

/**
 * Elementary Monte Carlo trial that exchanges volume between two phases.  Trial
 * consists of a volume increase in one phase (selected at random) and an equal
 * volume decrease in the other.  Used in Gibbs ensemble simulations.
 *
 * @author David Kofke
 */
public final class MCMoveVolumeExchange extends MCMove {
    
    private Phase firstPhase;
    private Phase secondPhase;
    private PhaseAction.Inflate inflate1;
    private PhaseAction.Inflate inflate2;
    private final double ROOT;
    private final IteratorDirective iteratorDirective = new IteratorDirective();
    private final PotentialCalculation.EnergySum energy = new PotentialCalculation.EnergySum();

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
    
    public void thisTrial() {
        double hOld = firstPhase.potential().calculate(iteratorDirective, energy.reset()).sum()
                    + secondPhase.potential().calculate(iteratorDirective, energy.reset()).sum();
        double v1Old = firstPhase.volume();
        double v2Old = secondPhase.volume();
        double vRatio = v1Old/v2Old * Math.exp(stepSize*(Simulation.random.nextDouble() - 0.5));
        double v2New = (v1Old + v2Old)/(1 + vRatio);
        double v1New = (v1Old + v2Old - v2New);
        double v1Scale = v1New/v1Old;
        double v2Scale = v2New/v2Old;
//        System.out.println(v1Scale + "  " + v2Scale + " " +stepSize + " "+nTrials+" "+nAccept+" "+adjustInterval+" "+frequency);
        inflate1.setScale(Math.pow(v1Scale,ROOT));
        inflate2.setScale(Math.pow(v2Scale,ROOT));
        inflate1.attempt();
        inflate2.attempt();
        double hNew = firstPhase.potential().calculate(iteratorDirective, energy.reset()).sum()
                    + secondPhase.potential().calculate(iteratorDirective, energy.reset()).sum();
        if(hNew >= Double.MAX_VALUE ||
             Math.exp(-(hNew-hOld)/parentIntegrator.temperature+
                       (firstPhase.moleculeCount()+1)*Math.log(v1Scale) +
                       (secondPhase.moleculeCount()+1)*Math.log(v2Scale))
                < Simulation.random.nextDouble()) 
            {  //reject
              inflate1.undo();
              inflate2.undo();
        }
        else nAccept++;
    }//end of thisTrial
}//end of MCMoveVolumeExchange