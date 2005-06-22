package etomica.data.meter;
import etomica.Phase;
import etomica.data.DataSourceScalar;
import etomica.units.Dimension;
import etomica.units.DimensionRatio;
	
/**
 * Class used to construct meters for intensive properties formed as the ratio of two extensive properties (e.g., energy/mole)
 */
public class MeterRatio extends DataSourceScalar {
	    
	protected DataSourceScalar nMeter;
	protected DataSourceScalar dMeter;
	    
	public MeterRatio(DataSourceScalar numerator, DataSourceScalar denominator) {
	    super();
	    setNumerator(numerator);
	    setDenominator(denominator);
	}
	    
	public Dimension getDimension() {
	    return new DimensionRatio(nMeter.getDimension(), dMeter.getDimension());
	}
	    
	public double getDataAsScalar(Phase p) {
        if (phase == null) throw new IllegalStateException("must call setPhase before using meter");
	    return nMeter.getDataAsScalar(p)/dMeter.getDataAsScalar(p);
	}
	
	public void setPhase(Phase[] p) {
	    super.setPhase(p);
	    nMeter.setPhase(p);
	    dMeter.setPhase(p);
	}
	    
	public void setNumerator(DataSourceScalar numerator) {nMeter = numerator;}
	public void setDenominator(DataSourceScalar denominator) {dMeter = denominator;}
	    
	    
}
