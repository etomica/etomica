package etomica;
import etomica.units.Dimension;

/**
 * Meter for evaluation of the soft-potential pressure in a phase.
 * Requires that temperature be set in order to calculation ideal-gas
 * contribution to pressure; default is to use zero temperature, which
 * causes this contribution to be omitted.
 *
 * @author David Kofke
 */
 
public class MeterPressure extends MeterScalar implements EtomicaElement {
    
    private IteratorDirective iteratorDirective;
    private final PotentialCalculationVirialSum virial;
    private final PotentialMaster potential;
    private double temperature;
    private final double rD;
    
    public MeterPressure() {
        this(Simulation.instance);
    }

    public MeterPressure(Simulation sim) {
        super(sim);
        setTemperature(temperature);
        rD = 1.0/(double)simulation().space.D();
        setLabel("Pressure");
        iteratorDirective = new IteratorDirective();
        iteratorDirective.includeLrc = true;
        potential = simulation().potentialMaster;
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
     * Returns Dimension.PRESSURE
     */
    public Dimension getDimension() {return Dimension.PRESSURE;}

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
    public double getDataAsScalar(Phase phase) {
    	virial.zeroSum();
        potential.calculate(phase, iteratorDirective, virial);
        return phase.getDensity()*temperature - virial.getSum()*rD/phase.boundary().volume();
    }
    
}//end of MeterPressure
