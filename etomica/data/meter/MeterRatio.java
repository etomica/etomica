package etomica.data.meter;
import etomica.Phase;
import etomica.units.Dimension;
import etomica.units.DimensionRatio;
	
/**
 * Class used to construct meters for intensive properties formed as the ratio of two extensive properties (e.g., energy/mole)
 */
public class MeterRatio extends MeterScalar {
	    
	protected MeterScalar nMeter;
	protected MeterScalar dMeter;
	    
	public MeterRatio(MeterScalar numerator, MeterScalar denominator) {
	    super();
	    setNumerator(numerator);
	    setDenominator(denominator);
	}
	    
	public Dimension getDimension() {
	    return new DimensionRatio(nMeter.getDimension(), dMeter.getDimension());
	}
	    
	public double getDataAsScalar(Phase p) {
	    return nMeter.getDataAsScalar(p)/dMeter.getDataAsScalar(p);
	}
	
	public void setPhase(Phase[] p) {
	    super.setPhase(p);
	    nMeter.setPhase(p);
	    dMeter.setPhase(p);
	}
	    
	public void setNumerator(MeterScalar numerator) {nMeter = numerator;}
	public void setDenominator(MeterScalar denominator) {dMeter = denominator;}
	    
	    
}
