package etomica.data;

import etomica.Constants;

/**
 * Accumulator for calculating ratio between two sums
 */
public class AccumulatorRatioAverage extends AccumulatorAverage {
    
    /**
     * Compute and return the ratio of the sums
     */
    public double[] ratio() {
        if(sum[0] == 0.0) setNaN(average);
        else {
	       	for(int i=nDataMinus1; i>=0; i--) {       		
	        	average[i] = sum[i]/sum[0];
	        }
        }
        return average;
    }
	
	/**
	 * Return the 67% confidence limits of the ratio based on variance in block averages.
	 */
    public double[] ratioError() {
    	if(count > 1) {
            // calculate actual average and error for each data set
            super.average();
            error();
            if (average[0] == 0.0) {
                setNaN(ratioError);
            }
            else {
                ratioError[0] = error[0] / average[0];
                ratioError[0] = ratioError[0] * ratioError[0];
                for(int i=nDataMinus1; i>=1; i--) {
                    if (average[i] == 0.0) {
                        ratioError[i] = Double.NaN;
                    }
                    else {
                        ratioError[i] = error[i]/average[i];
                        ratioError[i] = Math.sqrt(ratioError[i]*ratioError[i] + ratioError[0]) * average[i]/average[0];
                    }
                }
                ratioError[0] = 0.0;
            }
    	} else {
    		setNaN(ratioError);
    	}
        return ratioError;
    }
    
    /**
     * Compute and return the standard deviation of the ratio.
     */
     public double[] ratioStandardDeviation() {
     	if(count > 1) {
            // calculate actual average and standard deviation for each data set
            super.average();
            standardDeviation();
            if (average[0] == 0.0) {
                setNaN(ratioStandardDeviation);
            }
            else {
                ratioStandardDeviation[0] = standardDeviation[0] / average[0];
                ratioStandardDeviation[0] = ratioStandardDeviation[0] * ratioStandardDeviation[0];
                for(int i=nDataMinus1; i>=1; i--) {
                    if (average[i] == 0.0) {
                        ratioStandardDeviation[i] = Double.NaN;
                    }
                    else {
                        ratioStandardDeviation[i] = standardDeviation[i]/average[i];
                        ratioStandardDeviation[i] = Math.sqrt(ratioStandardDeviation[i]*ratioStandardDeviation[i] + ratioStandardDeviation[0]) * average[i]/average[0];
                    }
                }
                ratioStandardDeviation[0] = 0.0;
            }
        } else {
            setNaN(ratioStandardDeviation);
        }
     	return ratioStandardDeviation;
     }
    
    /**
     * Returns the value indicated by the argument.
     */
    public double[] getData(AccumulatorAverage.Type type) {
        if(type==RATIO) return ratio();
        if(type==RATIO_ERROR) return ratioError();
        if(type==RATIO_STANDARD_DEVIATION) return ratioStandardDeviation();
        return super.getData(type);
    }
    public static class Type extends AccumulatorAverage.Type {
        protected Type(String label) {super(label);}       
        public Constants.TypedConstant[] choices() {return VIRIAL_CHOICES;}
    }
    //XXX such an ugly hack!!!!
    protected static final AccumulatorAverage.Type[] VIRIAL_CHOICES = 
        new AccumulatorAverage.Type[] {
            CHOICES[0], CHOICES[1], CHOICES[2], CHOICES[3], CHOICES[4],
            new Type("Ratio"),
            new Type("Ratio error"), new Type("Ratio standard deviation")};
    public static final AccumulatorAverage.Type RATIO = VIRIAL_CHOICES[5];
    public static final AccumulatorAverage.Type RATIO_ERROR = VIRIAL_CHOICES[6];
    public static final AccumulatorAverage.Type RATIO_STANDARD_DEVIATION = VIRIAL_CHOICES[7];

    protected void setNData(int nData) {
        ratioStandardDeviation = redimension(nData, ratioError);
        ratioError = redimension(nData, ratioError);
    	 super.setNData(nData);
    }
     
    //need separate fields because ratio values are calculated from the non-ratio values.
    protected double[] ratioStandardDeviation, ratioError;
}
