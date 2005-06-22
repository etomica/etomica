package etomica.integrator.mcmove;

import etomica.AtomIterator;
import etomica.Phase;
import etomica.PotentialMaster;
import etomica.Simulation;
import etomica.Space;
import etomica.action.PhaseInflate;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.MCMove;

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

public final class MCMoveVolumeExchange extends MCMove {
    
    private final MeterPotentialEnergy energyMeter;
    private Phase firstPhase;
    private Phase secondPhase;
    private final PhaseInflate inflate1;
    private final PhaseInflate inflate2;
    private transient double uOld1, uOld2;
    private transient double uNew1 = Double.NaN;
    private transient double uNew2 = Double.NaN;
    private final double ROOT;
    private final AtomIteratorAllMolecules phase1AtomIterator;
    private final AtomIteratorAllMolecules phase2AtomIterator;
    
    private transient double hOld, v1Scale, v2Scale;

    public MCMoveVolumeExchange(PotentialMaster potentialMaster, Space space) {
        super(potentialMaster, 2);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        ROOT = 1.0/space.D();
        setStepSizeMax(Double.MAX_VALUE);
        setStepSizeMin(Double.MIN_VALUE);
        setStepSize(0.3);
        phase1AtomIterator = new AtomIteratorAllMolecules();
        phase2AtomIterator = new AtomIteratorAllMolecules();
        energyMeter.setIncludeLrc(false);
        inflate1 = new PhaseInflate(space);
        inflate2 = new PhaseInflate(space);
    }
    
    public void setPhase(Phase[] p) {
        super.setPhase(p);
        firstPhase = p[0];
        secondPhase = p[1];
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
        double vRatio = v1Old/v2Old * Math.exp(stepSize*(Simulation.random.nextDouble() - 0.5));
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
    
    public double lnTrialRatio() {
        return (firstPhase.moleculeCount()+1)*Math.log(v1Scale) +
                + (secondPhase.moleculeCount()+1)*Math.log(v2Scale);
    }
        
    public double lnProbabilityRatio() {
        energyMeter.setPhase(firstPhase);
        uNew1 = energyMeter.getDataAsScalar();
        energyMeter.setPhase(secondPhase);
        uNew2 = energyMeter.getDataAsScalar();
        double hNew = uNew1 + uNew2;
        return -(hNew - hOld)/temperature;
    }
    
    public void acceptNotify() {  /* do nothing */}
    
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
            return AtomIterator.NULL;
        }
    }

}//end of MCMoveVolumeExchange