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
    private AtomIterator phase1AtomIterator;
    private AtomIterator phase2AtomIterator;
//    private AtomIteratorCompound affectedAtomIterator;
    
    private transient double hOld, v1Scale, v2Scale;

    public MCMoveVolumeExchange(IntegratorMC parent) {
        super(parent);
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
        if(p == null || p.length == 0) return;
        super.setPhase(p);
        firstPhase = p[0];
        if(p.length < 2) return;
        secondPhase = p[1];
        if(firstPhase == null && secondPhase == null) return;
        inflate1 = new PhaseAction.Inflate(firstPhase);
        inflate2 = new PhaseAction.Inflate(secondPhase);
        phase1AtomIterator = firstPhase.makeMoleculeIterator();
        phase2AtomIterator = secondPhase.makeMoleculeIterator();
//        affectedAtomIterator 
//            = new AtomIteratorCompound(new AtomIterator[] {phase1AtomIterator, phase2AtomIterator});
    }
    
    public boolean doTrial() {
        hOld = potential.set(firstPhase).calculate(iteratorDirective, energy.reset()).sum()
                    + potential.set(secondPhase).calculate(iteratorDirective, energy.reset()).sum();
        double v1Old = firstPhase.volume();
        double v2Old = secondPhase.volume();
        double vRatio = v1Old/v2Old * Math.exp(stepSize*(Simulation.random.nextDouble() - 0.5));
        double v2New = (v1Old + v2Old)/(1 + vRatio);
        double v1New = (v1Old + v2Old - v2New);
        v1Scale = v1New/v1Old;
        v2Scale = v2New/v2Old;
        inflate1.setScale(Math.pow(v1Scale,ROOT));
        inflate2.setScale(Math.pow(v2Scale,ROOT));
        inflate1.attempt();
        inflate2.attempt();
        return true;
    }//end of doTrial
    
    public double lnTrialRatio() {
        return (firstPhase.moleculeCount()+1)*Math.log(v1Scale) +
                + (secondPhase.moleculeCount()+1)*Math.log(v2Scale);
    }
        
    public double lnProbabilityRatio() {
        double hNew = potential.set(firstPhase).calculate(iteratorDirective, energy.reset()).sum()
                    + potential.set(secondPhase).calculate(iteratorDirective, energy.reset()).sum();
        return -(hNew - hOld)/parentIntegrator.temperature;
    }
    
    public void acceptNotify() {  /* do nothing */}
    
    public void rejectNotify() {
        inflate1.undo();
        inflate2.undo();
    }

    public final AtomIterator affectedAtoms(Phase phase) {
        if(this.firstPhase == phase) {
            phase1AtomIterator.reset();
            return phase1AtomIterator;
        } else if(this.secondPhase == phase) {
            phase2AtomIterator.reset();
            return phase2AtomIterator;
        } else {
            return AtomIterator.NULL;
        }
    }

}//end of MCMoveVolumeExchange