package etomica.data.meter;
import etomica.AtomSet;
import etomica.EtomicaInfo;
import etomica.IteratorDirective;
import etomica.Meter;
import etomica.Phase;
import etomica.PotentialMaster;
import etomica.data.DataSourceScalar;
import etomica.potential.PotentialCalculationEnergySum;
import etomica.units.Dimension;

/**
 * Meter for evaluation of the potential energy in a phase.
 * Includes several related methods for computing the potential energy of a single
 * atom or molecule with all neighboring atoms
 *
 * @author David Kofke
 */
 
public class MeterPotentialEnergy extends DataSourceScalar implements Meter {
    
    public MeterPotentialEnergy(PotentialMaster potentialMaster) {
        super("Potential Energy",Dimension.ENERGY);
        iteratorDirective.includeLrc = true;
        potential = potentialMaster;
        iteratorDirective.setDirection(null); // so that "both" will work
    }
      
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Total intermolecular potential energy in a phase");
        return info;
    }

    /**
     * Sets flag indicating whether calculated energy should include
     * long-range correction for potential truncation (true) or not (false).
     */
    public void setIncludeLrc(boolean b) {
    	iteratorDirective.includeLrc = b;
    }
    /**
     * Indicates whether calculated energy should include
     * long-range correction for potential truncation (true) or not (false).
     */
    public boolean isIncludeLrc() {
    	return iteratorDirective.includeLrc;
    }

    public void setTarget(AtomSet atoms) {
    	iteratorDirective.setTargetAtoms(atoms);
    }
    
   /**
    * Computes total potential energy for phase.
    * Currently, does not include long-range correction to truncation of energy
    */
    public double getDataAsScalar() {
        if (phase == null) throw new IllegalStateException("must call setPhase before using meter");
    	energy.zeroSum();
        potential.calculate(phase, iteratorDirective, energy);
        return energy.getSum();
    }
    /**
     * @return Returns the phase.
     */
    public Phase getPhase() {
        return phase;
    }
    /**
     * @param phase The phase to set.
     */
    public void setPhase(Phase phase) {
        this.phase = phase;
    }

    private Phase phase;
    private final IteratorDirective iteratorDirective = new IteratorDirective();
    private final PotentialCalculationEnergySum energy = new PotentialCalculationEnergySum();
    private final PotentialMaster potential;
    
}//end of MeterPotentialEnergy
