package etomica.data;

import etomica.data.types.DataArithmetic;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.simulation.Simulation;
import etomica.units.Null;
import etomica.util.EnumeratedType;
import etomica.util.Function;

/**
 * Accumulator for calculating ratio between two sums
 */
public class AccumulatorRatioAverage extends AccumulatorAverage {
    
    public AccumulatorRatioAverage(Simulation sim) {
        this(sim.getDefaults().blockSize);
    }
    
    public AccumulatorRatioAverage(int blockSize) {
        super(blockSize);
    }
    
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
            ratioError.DE(average);
            ratioError.TE(ratioError);
            ratioError.PE(errorRatio0);
            ratioError.map(Function.Sqrt.INSTANCE);
            ratioError.TE(ratio);
            ratioError.map(Function.Abs.INSTANCE);

            double stdevRatio0 = ((DataDoubleArray)standardDeviation).getData()[0]/average0;
            ratioStandardDeviation.E((Data)standardDeviation);
            ratioStandardDeviation.DE(average);
            ratioStandardDeviation.TE(ratioStandardDeviation);
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
    
    public DataInfo processDataInfo(DataInfo dataInfo) {
        DataFactory factory = dataInfo.getDataFactory();
        
        sum = (DataArithmetic)factory.makeData(dataInfo.getLabel()+"(blkAvg sum)", dataInfo.getDimension());

        ratio = (DataArithmetic)factory.makeData("Ratio", Null.DIMENSION);
        ratioError = (DataArithmetic)factory.makeData("Ratio error", Null.DIMENSION);
        ratioStandardDeviation = (DataArithmetic)factory.makeData("Ratio stddev", Null.DIMENSION);
        super.processDataInfo(dataInfo);
        Data[] dataGroups = new Data[dataGroup.getNData()+3];
        int i;
        for (i=0; i<dataGroup.getNData(); i++) {
            dataGroups[i] = dataGroup.getData(i);
        }
        dataGroups[i++] = (Data)ratio;
        dataGroups[i++] = (Data)ratioError;
        dataGroups[i++] = (Data)ratioStandardDeviation;
        dataGroup = new DataGroup("Group", dataGroups);
        return dataGroup.getDataInfo();
    }
    
    public static class StatType extends AccumulatorAverage.StatType {
        protected StatType(String label, int index) {super(label,index);}       
        public static final StatType RATIO = new StatType("Ratio",6);
        public static final StatType RATIO_ERROR = new StatType("Ratio error",7);
        public static final StatType RATIO_STANDARD_DEVIATION = new StatType("Ratio standard deviation",8);
        public static AccumulatorAverage.StatType[] choices() {
            AccumulatorAverage.StatType[] choices = AccumulatorAverage.StatType.choices();
            return new AccumulatorAverage.StatType[] {
                choices[0], choices[1], choices[2], choices[3], choices[4], choices[5],
                RATIO, RATIO_ERROR, RATIO_STANDARD_DEVIATION};
        }
    }

    //need separate fields because ratio values are calculated from the non-ratio values.
    protected DataArithmetic ratio, ratioStandardDeviation, ratioError;
}
