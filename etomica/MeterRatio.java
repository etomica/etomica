package etomica;
import etomica.units.*;
	
/**
 * Class used to construct meters for intensive properties formed as the ratio of two extensive properties (e.g., energy/mole)
 */
public class MeterRatio extends Meter {
	    
	protected Meter nMeter;
	protected Meter dMeter;
	    
	public MeterRatio(Meter numerator, Meter denominator) {
	    this(Simulation.instance, numerator, denominator);
	}
	public MeterRatio(Simulation sim, Meter numerator, Meter denominator) {
	    super(sim);
	    setNumerator(numerator);
	    setDenominator(denominator);
	}
	    
    /**
     * Declaration whether this meter uses the boundary object of phase when making its measurements
     */
    public final boolean usesPhaseBoundary() {
        return nMeter.usesPhaseBoundary() || dMeter.usesPhaseBoundary();
    }
    /**
     * Declaration whether this meter uses the iteratorFactory of phase when making its measurements
     */
    public final boolean usesPhaseIteratorFactory() {
        return nMeter.usesPhaseIteratorFactory() || dMeter.usesPhaseIteratorFactory();
    }


	public Dimension getDimension() {
	    return new DimensionRatio(nMeter.getDimension(), dMeter.getDimension());
	}
	    
	public double currentValue() {
	    return nMeter.currentValue()/dMeter.currentValue();
	}
	    
	public void setNumerator(Meter numerator) {nMeter = numerator;}
	public void setDenominator(Meter denominator) {dMeter = denominator;}
	    
	    
}
