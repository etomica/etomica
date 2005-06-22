package etomica.data;
import etomica.DataInfo;
import etomica.units.DimensionRatio;
	
/**
 * Class used to construct meters for intensive properties formed as the ratio of two extensive properties (e.g., energy/mole)
 */
public class DataSourceRatio extends DataSourceScalar {
	    
	protected DataSourceScalar nMeter;
	protected DataSourceScalar dMeter;
	    
	public DataSourceRatio(DataSourceScalar numerator, DataSourceScalar denominator) {
	    super(new DataInfo("Ratio "+numerator.getDataInfo().getLabel()+"/"+denominator.getDataInfo().getLabel(),
                new DimensionRatio(numerator.getDataInfo().getDimension(),denominator.getDataInfo().getDimension())));
	    setNumerator(numerator);
	    setDenominator(denominator);
	}
	    
	public double getDataAsScalar() {
	    return nMeter.getDataAsScalar()/dMeter.getDataAsScalar();
	}
	
	public void setNumerator(DataSourceScalar numerator) {nMeter = numerator;}
	public void setDenominator(DataSourceScalar denominator) {dMeter = denominator;}
}
