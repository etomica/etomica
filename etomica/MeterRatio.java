package etomica;
import etomica.units.*;
	
/**
 * Class used to construct meters for intensive properties formed as the ratio of two extensive properties (e.g., energy/mole)
 */
public class MeterRatio extends MeterScalar {
	    
	protected MeterScalar nMeter;
	protected MeterScalar dMeter;
	    
	public MeterRatio(MeterScalar numerator, MeterScalar denominator) {
	    this(Simulation.instance, numerator, denominator);
	}
	public MeterRatio(Simulation sim, MeterScalar numerator, MeterScalar denominator) {
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
	    
	public void setNumerator(MeterScalar numerator) {nMeter = numerator;}
	public void setDenominator(MeterScalar denominator) {dMeter = denominator;}
	    
	    
}
