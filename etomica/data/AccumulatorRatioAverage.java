package etomica.data;

import etomica.Constants;
import etomica.Data;
import etomica.data.types.DataArithmetic;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.utility.Function;

/**
 * Accumulator for calculating ratio between two sums
 */
public class AccumulatorRatioAverage extends AccumulatorAverage {
    
    public Data getData() {
        if (mostRecent == null) return null;
        if(count > 0) {
            super.getData();
            double average0 = ((DataDoubleArray)average).getData()[0];
            if (average0 == 0) {
                ratio.E(Double.NaN);
                ratioError.E(Double.NaN);
                ratioStandardDeviation.E(Double.NaN);
                return dataGroup;
            }

            ratio.E((Data)sum);
            ratio.TE(1/((DataDoubleArray)sum).getData()[0]);

            double errorRatio0 = ((DataDoubleArray)error).getData()[0]/average0;
            errorRatio0 *= errorRatio0;
            ratioError.E((Data)error);
            ratioError.TE(error);
            ratioError.PE(errorRatio0);
            ratioError.map(Function.Sqrt.INSTANCE);
            ratioError.TE(ratio);
            ratioError.map(Function.Abs.INSTANCE);

            double stdevRatio0 = ((DataDoubleArray)standardDeviation).getData()[0]/average0;
            ratioStandardDeviation.E((Data)standardDeviation);
            ratioStandardDeviation.TE(standardDeviation);
            ratioStandardDeviation.PE(stdevRatio0);
            ratioStandardDeviation.map(Function.Sqrt.INSTANCE);
            ratioStandardDeviation.TE(ratio);
            ratioStandardDeviation.map(Function.Abs.INSTANCE);
        }
        return dataGroup;
    }

    public void reset() {
        super.reset();
        ratio.E(Double.NaN);
        ratioError.E(Double.NaN);
        ratioStandardDeviation.E(Double.NaN);
    }
    
    protected void initialize(Data value) {
        ratio = (DataArithmetic)value.makeCopy();
        ratioError = (DataArithmetic)value.makeCopy();
        ratioStandardDeviation = (DataArithmetic)value.makeCopy();
        super.initialize(value);
        Data[] dataGroups = new Data[dataGroup.getNData()+2];
        int i;
        for (i=0; i<dataGroup.getNData(); i++) {
            dataGroups[i] = dataGroup.getData(i);
        }
        dataGroups[i++] = (Data)ratio;
        dataGroups[i++] = (Data)ratioError;
        dataGroups[i++] = (Data)ratioStandardDeviation;
        dataGroup = new DataGroup(value.getDataInfo(),dataGroups);
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

    //need separate fields because ratio values are calculated from the non-ratio values.
    protected DataArithmetic ratio, ratioStandardDeviation, ratioError;
}
