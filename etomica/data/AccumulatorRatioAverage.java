package etomica.data;

import etomica.Constants;

/**
 * Accumulator for calculating ratio between two sums
 */
public class AccumulatorRatioAverage extends AccumulatorAverage {
    
    public double[] getData() {
        int currentBlockCount = blockSize - blockCountDown;
        if(count+currentBlockCount == 0) {
            setNaN(data);
        } else {
            super.getData();
            int numBaseStats = super.numStats();
            data[(numBaseStats+0)*nData+0] = 1; // average
            data[(numBaseStats+1)*nData+0] = 0; // error
            data[(numBaseStats+2)*nData+0] = 0; // std dev
            double errorRatio0 = error[0]/average[0];
            errorRatio0 *= errorRatio0;
            double stdevRatio0 = standardDeviation[0]/average[0];
            stdevRatio0 *= stdevRatio0;
            for(int i=1; i<nData; i++) {
                double errorRatio = Double.NaN;
                double stdevRatio = Double.NaN;
                if (average[0] != 0.0) {
                    errorRatio = error[i]/average[i];
                    errorRatio = Math.sqrt(errorRatio*errorRatio + errorRatio0) * Math.abs(average[i]/average[0]);
                    stdevRatio = standardDeviation[i]/average[i];
                    stdevRatio = Math.sqrt(stdevRatio*stdevRatio + stdevRatio0) * Math.abs(average[i]/average[0]);
                }
                data[(numBaseStats+0)*nData+i] = sum[i]/sum[0];
                data[(numBaseStats+1)*nData+i] = ratioError[i] = errorRatio;
                data[(numBaseStats+2)*nData+i] = ratioStandardDeviation[i] = stdevRatio;
            }
        }
        return data;
    }
    
    public int numStats() {
        return 3+super.numStats();
    }
     
    public static class Type extends AccumulatorAverage.Type {
        protected Type(String label, int index) {super(label,index);}       
        public Constants.TypedConstant[] choices() {return VIRIAL_CHOICES;}
    }
    //XXX such an ugly hack!!!!
    protected static final AccumulatorAverage.Type[] VIRIAL_CHOICES = 
        new AccumulatorAverage.Type[] {
            CHOICES[0], CHOICES[1], CHOICES[2], CHOICES[3], CHOICES[4],
            new Type("Ratio",5),
            new Type("Ratio error",6), new Type("Ratio standard deviation",7)};
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
