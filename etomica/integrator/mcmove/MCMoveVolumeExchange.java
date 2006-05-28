package etomica.integrator.mcmove;

import etomica.action.PhaseInflate;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.atom.iterator.AtomIteratorNull;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.IntegratorPhase;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;

/**
 * Elementary Monte Carlo trial that exchanges volume between two phases.  Trial
 * consists of a volume increase in one phase (selected at random) and an equal
 * volume decrease in the other.  Used in Gibbs ensemble simulations.
 *
 * @author David Kofke
 */
 
 /* History of changes
  * 7/9/02 added energyChange() method.
  */

public final class MCMoveVolumeExchange extends MCMoveStep {
    
    private final MeterPotentialEnergy energyMeter;
    private final Phase firstPhase;
    private final Phase secondPhase;
    private final IntegratorPhase integrator1;
    private final IntegratorPhase integrator2;
    private final PhaseInflate inflate1;
    private final PhaseInflate inflate2;
    private transient double uOld1, uOld2;
    private transient double uNew1 = Double.NaN;
    private transient double uNew2 = Double.NaN;
    private final double ROOT;
    private final AtomIteratorAllMolecules phase1AtomIterator;
    private final AtomIteratorAllMolecules phase2AtomIterator;
    
    private transient double hOld, v1Scale, v2Scale;

    public MCMoveVolumeExchange(PotentialMaster potentialMaster,
            IntegratorPhase integrator1, IntegratorPhase integrator2) {
        super(potentialMaster, new MCMoveStepTracker());
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        Space space = potentialMaster.getSpace();
        ROOT = 1.0/space.D();
        setStepSizeMax(Double.MAX_VALUE);
        setStepSizeMin(Double.MIN_VALUE);
        setStepSize(0.3);
        phase1AtomIterator = new AtomIteratorAllMolecules();
        phase2AtomIterator = new AtomIteratorAllMolecules();
        energyMeter.setIncludeLrc(false);
        inflate1 = new PhaseInflate(space);
        inflate2 = new PhaseInflate(space);
        this.integrator1 = integrator1;
        this.integrator2 = integrator2;
        firstPhase = integrator1.getPhase();
        secondPhase = integrator2.getPhase();
        inflate1.setPhase(firstPhase);
        inflate2.setPhase(secondPhase);
        phase1AtomIterator.setPhase(firstPhase);
        phase2AtomIterator.setPhase(secondPhase);
    }
    
    public boolean doTrial() {
        energyMeter.setPhase(firstPhase);
        uOld1 = energyMeter.getDataAsScalar();
        energyMeter.setPhase(secondPhase);
        uOld2 = energyMeter.getDataAsScalar();
        hOld = uOld1 + uOld2;
        double v1Old = firstPhase.volume();
        double v2Old = secondPhase.volume();
        double step = stepSize * (Simulation.random.nextDouble() - 0.5); 
        double vRatio = v1Old/v2Old * Math.exp(step);
        double v2New = (v1Old + v2Old)/(1 + vRatio);
        double v1New = (v1Old + v2Old - v2New);
        v1Scale = v1New/v1Old;
        v2Scale = v2New/v2Old;
        inflate1.setScale(Math.pow(v1Scale,ROOT));
        inflate2.setScale(Math.pow(v2Scale,ROOT));
        inflate1.actionPerformed();
        inflate2.actionPerformed();
        uNew1 = uNew2 = Double.NaN;
        return true;
    }//end of doTrial
    
    public double getA() {
        energyMeter.setPhase(firstPhase);
        uNew1 = energyMeter.getDataAsScalar();
        energyMeter.setPhase(secondPhase);
        uNew2 = energyMeter.getDataAsScalar();
        double hNew = uNew1 + uNew2;
        double B = -(hNew - hOld);
        // assume both integrators have the same temperature
        double T = integrator1.getTemperature();
        return Math.exp(B/T) * Math.pow(v1Scale,(firstPhase.moleculeCount()+1))
                * Math.pow(v2Scale,(secondPhase.moleculeCount()+1));
    }
        
    public double getB() {
        //IntegratorManagerMC only calls getA since it doesn't have a temperature
        throw new IllegalStateException("You shouldn't be calling this method");
    }
    
    public void acceptNotify() {
        try {
            integrator1.reset();
            integrator2.reset();
        } catch(ConfigurationOverlapException e) {
            throw new RuntimeException(e);
        }
    }
    
    public void rejectNotify() {
        inflate1.undo();
        inflate2.undo();
    }

    public double energyChange(Phase phase) {
        if(this.firstPhase == phase) return uNew1 - uOld1;
        else if(this.secondPhase == phase) return uNew2 - uOld2;
        else return 0.0;
    }
    
    public final AtomIterator affectedAtoms(Phase phase) {
        if(this.firstPhase == phase) {
            return phase1AtomIterator;
        } else if(this.secondPhase == phase) {
            return phase2AtomIterator;
        } else {
            return AtomIteratorNull.INSTANCE;
        }
    }

}//end of MCMoveVolumeExchange