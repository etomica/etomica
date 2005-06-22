package etomica.data.meter;
import etomica.DataInfo;
import etomica.EtomicaInfo;
import etomica.IteratorDirective;
import etomica.Meter;
import etomica.Phase;
import etomica.PotentialMaster;
import etomica.Space;
import etomica.data.DataSourceScalar;
import etomica.potential.PotentialCalculationVirialSum;
import etomica.units.Dimension;

/**
 * Meter for evaluation of the soft-potential pressure in a phase.
 * Requires that temperature be set in order to calculation ideal-gas
 * contribution to pressure; default is to use zero temperature, which
 * causes this contribution to be omitted.
 *
 * @author David Kofke
 */
 
public class MeterPressure extends DataSourceScalar implements Meter {
    
    private IteratorDirective iteratorDirective;
    private final PotentialCalculationVirialSum virial;
    private final PotentialMaster potential;
    private double temperature;
    private final double rD;
    
    public MeterPressure(PotentialMaster potentialMaster, Space space) {
    	super(new DataInfo("Pressure",Dimension.pressure(space.D)));
        setTemperature(temperature);
        rD = 1.0/space.D();
        iteratorDirective = new IteratorDirective();
        iteratorDirective.includeLrc = true;
        potential = potentialMaster;
        virial = new PotentialCalculationVirialSum();
    }
      
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Total pressure in a phase (requires soft-potential model)");
        return info;
    }

    /**
     * Returns the temperature used to compute the ideal-gas contribution to
     * the pressure.
     * @return
     */
	public double getTemperature() {
		return temperature;
	}
	/**
	 * Sets the temperature used to compute the ideal-gas contribution to the pressure.
	 * @param temperature
	 */
	public void setTemperature(double temperature) {
		this.temperature = temperature;
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

	 /**
	  * Computes total pressure in phase by summing virial over all pairs, and adding
	  * ideal-gas contribution.
	  * Currently, does not include long-range correction to truncation of energy.
	  */
    public double getDataAsScalar() {
        if (phase == null) throw new IllegalStateException("must call setPhase before using meter");
    	virial.zeroSum();
        potential.calculate(phase, iteratorDirective, virial);
        return phase.getDensity()*temperature - virial.getSum()*rD/phase.boundary().volume();
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
}//end of MeterPressure
