package etomica.data;

import etomica.data.types.DataGroup;
import etomica.data.types.DataGroup.DataInfoGroupFactory;
import etomica.util.Arrays;
import etomica.util.Function;

/**
 * Accumulator for calculating ratio between two sums
 */
public class AccumulatorRatioAverageCovariance extends AccumulatorAverageCovariance {
    
    private static final long serialVersionUID = 1L;
    public AccumulatorRatioAverageCovariance() {
        this(1000);
    }
    
    public AccumulatorRatioAverageCovariance(long blockSize) {
        super(blockSize);
        ratioTag = new DataTag();
        ratioStandardDeviationTag = new DataTag();
        ratioErrorTag = new DataTag();
    }
    
    public DataTag getTag(StatType statType) {
        if (statType == RATIO) {
            return ratioTag;
        }
        if (statType == RATIO_STANDARD_DEVIATION) {
            return ratioStandardDeviationTag;
        }
        if (statType == RATIO_ERROR) {
            return ratioErrorTag;
        }
        return super.getTag(statType);
    }
    
    public IData getData() {
        if (sum == null) return null;
        if (count > 0) {
            super.getData();

            double average0 = average.getValue(0);
            if (average0 == 0) {
                ratio.E(Double.NaN);
                ratioError.E(Double.NaN);
                ratioStandardDeviation.E(Double.NaN);
                return dataGroup;
            }

            ratio.E(sum);
            ratio.TE(1/sum.getValue(0));


            work.E(sum);
            work.TE(1.0 / (count*blockSize));
            ratioError.E(error);
            ratioError.DE(work);
            ratioError.TE(ratioError);
            ratioError.PE(ratioError.getValue(0));
            otherErrorContribution = Math.sqrt(Math.abs(ratioError.getValue(1)))*Math.abs(ratio.getValue(1));
            double covarianceContribution = -2*blockCovariance.getValue(1)/(average.getValue(0)*average.getValue(1)*(count-1));
            covarianceErrorContribution = Math.sqrt(Math.abs(covarianceContribution))*Math.abs(ratio.getValue(1));
            ratioError.PE(covarianceContribution);
            ratioError.TE(ratio);
            ratioError.TE(ratio);
            ratioError.map(negativeChop);
            ratioError.map(Function.Sqrt.INSTANCE);
        }
        long nTotalData = count*blockSize + (blockSize-blockCountDown);
        if (nTotalData > 0) {
            // now use *all* of the data
            double average0 = average.getValue(0);
            ratio.E(average);
            ratio.TE(1/average0);

            double stdevRatio0 = standardDeviation.getValue(0)/average0;
            ratioStandardDeviation.E(standardDeviation);
            ratioStandardDeviation.DE(average);
            ratioStandardDeviation.TE(ratioStandardDeviation);
            ratioStandardDeviation.PE(stdevRatio0);
            ratioStandardDeviation.TE(ratio);
            ratioStandardDeviation.TE(ratio);
            ratioStandardDeviation.map(Function.Sqrt.INSTANCE);
        }

        return dataGroup;
    }

    public void reset() {
        if (sum == null || ratio == null) {
            //no data has been added yet, so nothing to reset
            return;
        }
        super.reset();
        ratio.E(Double.NaN);
        ratioError.E(Double.NaN);
        ratioStandardDeviation.E(Double.NaN);
    }
    
    public IEtomicaDataInfo processDataInfo(IEtomicaDataInfo incomingDataInfo) {
        super.processDataInfo(incomingDataInfo);

        ratio = incomingDataInfo.makeData();
        ratioError = incomingDataInfo.makeData();
        ratioStandardDeviation = incomingDataInfo.makeData();

        IData[] dataGroups = new IData[dataGroup.getNData()+3];
        int i;
        for (i=0; i<dataGroup.getNData(); i++) {
            dataGroups[i] = dataGroup.getData(i);
        }
        dataGroups[i++] = ratio;
        dataGroups[i++] = ratioError;
        dataGroups[i++] = ratioStandardDeviation;
        dataGroup = new DataGroup(dataGroups);
        reset();

        DataInfoGroupFactory groupFactory = (DataInfoGroupFactory)dataInfo.getFactory();
        IDataInfo[] subDataInfo = groupFactory.getSubDataInfo();

        IEtomicaDataInfoFactory factory = incomingDataInfo.getFactory();
        String incomingLabel = incomingDataInfo.getLabel();
        factory.setLabel(incomingLabel+" ratio");
        IEtomicaDataInfo ratioInfo = factory.makeDataInfo();
        ratioInfo.addTag(ratioTag);
        factory.setLabel(incomingLabel+" ratio error");
        IEtomicaDataInfo ratioErrorInfo = factory.makeDataInfo();
        ratioErrorInfo.addTag(ratioErrorTag);
        factory.setLabel(incomingLabel+" ratio");
        IEtomicaDataInfo ratioStandardDeviationInfo = factory.makeDataInfo();
        ratioStandardDeviationInfo.addTag(ratioStandardDeviationTag);
        
        subDataInfo = (IDataInfo[])Arrays.addObject(subDataInfo, ratioInfo);
        subDataInfo = (IDataInfo[])Arrays.addObject(subDataInfo, ratioErrorInfo);
        subDataInfo = (IDataInfo[])Arrays.addObject(subDataInfo, ratioStandardDeviationInfo);
        groupFactory.setSubDataInfo(subDataInfo);

        dataInfo = groupFactory.makeDataInfo();
        return dataInfo;
    }
    
    public static final StatType RATIO = new StatType("Ratio",7);
    public static final StatType RATIO_ERROR = new StatType("Ratio error",8);
    public static final StatType RATIO_STANDARD_DEVIATION = new StatType("Ratio standard deviation",9);
    public static AccumulatorAverage.StatType[] choices() {
        StatType[] choices = AccumulatorAverageCovariance.statChoices();
        return new AccumulatorAverage.StatType[] {
            choices[0], choices[1], choices[2], choices[3], choices[4], choices[5], choices[6],
            RATIO, RATIO_ERROR, RATIO_STANDARD_DEVIATION};
    }

    //need separate fields because ratio values are calculated from the non-ratio values.
    protected IData ratio, ratioStandardDeviation, ratioError;
    private final DataTag ratioTag, ratioStandardDeviationTag, ratioErrorTag;
    public double covarianceErrorContribution, otherErrorContribution;
}
