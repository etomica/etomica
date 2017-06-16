/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.data.types.DataGroup;
import etomica.data.types.DataGroup.DataInfoGroupFactory;
import etomica.util.Arrays;
import etomica.math.function.Function;

/**
 * Accumulator for calculating ratio between sums.  The ratios calculated are
 * the ratios of each data element to the first one, so v0/v0, v1/v0, v2/v0, etc.
 */
public class AccumulatorRatioAverageCovariance extends AccumulatorAverageCovariance {
    
    private static final long serialVersionUID = 1L;
    public AccumulatorRatioAverageCovariance() {
        this(1000);
    }
    
    public AccumulatorRatioAverageCovariance(long blockSize) {
        super(blockSize);
        ratioTag = new DataTag();
        ratioErrorTag = new DataTag();
    }
    
    public DataTag getTag(StatType statType) {
        if (statType == RATIO) {
            return ratioTag;
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
                return dataGroup;
            }

            ratio.E(sum);
            ratio.TE(1/sum.getValue(0));


            work.E(sum);
            work.TE(1.0 / (count*blockSize));
            ratioError.E(error);
            if (count > 1) {
                ratioError.DE(work);
                ratioError.TE(ratioError);
                ratioError.PE(ratioError.getValue(0));
                if (average.getValue(0)*average.getValue(1) > 0) {
                    double covarianceContribution = -2*blockCovariance.getValue(1)/(average.getValue(0)*average.getValue(1)*(count-1));
                    ratioError.PE(covarianceContribution);
                }
                ratioError.TE(ratio);
                ratioError.TE(ratio);
                ratioError.map(negativeChop);
                ratioError.map(Function.Sqrt.INSTANCE);
            }
        }
        long nTotalData = count*blockSize + (blockSize-blockCountDown);
        if (nTotalData > 0 && !doStrictBlockData) {
            if (count == 0) super.getData();
            // now use *all* of the data
            double average0 = average.getValue(0);
            ratio.E(average);
            ratio.TE(1/average0);
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
    }
    
    public IEtomicaDataInfo processDataInfo(IEtomicaDataInfo incomingDataInfo) {
        super.processDataInfo(incomingDataInfo);

        ratio = incomingDataInfo.makeData();
        ratioError = incomingDataInfo.makeData();

        IData[] dataGroups = new IData[dataGroup.getNData()+2];
        int i;
        for (i=0; i<dataGroup.getNData(); i++) {
            dataGroups[i] = dataGroup.getData(i);
        }
        dataGroups[i++] = ratio;
        dataGroups[i++] = ratioError;
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
        
        subDataInfo = (IDataInfo[])Arrays.addObject(subDataInfo, ratioInfo);
        subDataInfo = (IDataInfo[])Arrays.addObject(subDataInfo, ratioErrorInfo);
        groupFactory.setSubDataInfo(subDataInfo);

        dataInfo = groupFactory.makeDataInfo();
        return dataInfo;
    }
    
    public static final StatType RATIO = new StatType("Ratio",7);
    public static final StatType RATIO_ERROR = new StatType("Ratio error",8);
    public static AccumulatorAverage.StatType[] choices() {
        StatType[] choices = AccumulatorAverageCovariance.statChoices();
        return new AccumulatorAverage.StatType[] {
            choices[0], choices[1], choices[2], choices[3], choices[4], choices[5], choices[6],
            RATIO, RATIO_ERROR};
    }

    //need separate fields because ratio values are calculated from the non-ratio values.
    protected IData ratio, ratioError;
    private final DataTag ratioTag, ratioErrorTag;
}
