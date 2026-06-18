/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.data.types.DataGroup.DataInfoGroupFactory;
import etomica.units.dimensions.Null;
import etomica.util.Arrays;

/**
 * Accumulator for calculating ratio between sums.  This class calculates all
 * possible ratios vi/vj (as opposed to AccumulatorRatioAverageCovariance,
 * which calculates vi/v0.  
 */
public class AccumulatorRatioAverageCovarianceFull extends AccumulatorAverageCovariance {

    public AccumulatorRatioAverageCovarianceFull() {
        this(1000);
    }

    public AccumulatorRatioAverageCovarianceFull(long blockSize) {
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
        if (average == null) return null;
        // let superclass calculate its quantities
        super.getData();
        // now, calculate ratios
        int n = average.getLength();
        double[] r = ratio.getData();
        double[] re = ratioError.getData();
        if (count > 1) {
            // calculate uncertainties on ratios
            for (int i=0; i<n-1; i++) {
                re[i*n+i] = 0;
                double avgi = average.getValue(i);
                double rei = error.getValue(i)/(avgi);
                rei *= rei;
                for (int j=i+1; j<n; j++) {
                    double avgj = average.getValue(j);
                    double rej = error.getValue(j)/(avgj);
                    rej *= rej;
                    int k = i*n+j;
                    double rk = avgi/avgj;
                    double covarianceContribution = -2*blockCovariance.getValue(i*n+j)/(avgi*avgj*(count-1));
                    double rek = rk*rk * (rei + rej + covarianceContribution);
                    if (Double.isNaN(rek)) {
                        re[j*n+i] = Double.NaN;
                        re[k] = Double.NaN;
                        continue;
                    }
                    rek = (rek >= 0) ? rek : 0;
                    re[k] = Math.sqrt(rek);
                    
                    re[j*n+i] = re[k] / (rk*rk);
                }
            }
        }
        else {
        	ratioError.E(Double.NaN);
        }
        if (count > 0 || blockCountDown < blockSize) {
            // now use *all* of the data to calculate ratios
            for (int i=0; i<n; i++) {
                r[i*n+i] = 1;
                for (int j=i+1; j<n; j++) {
                    r[i*n+j] = average.getValue(i)/average.getValue(j);
                }
            }
        }

        return dataGroup;
    }

    public void reset() {
        if (average == null || ratio == null) {
            //no data has been added yet, so nothing to reset
            return;
        }
        super.reset();
        ratio.E(Double.NaN);
        ratioError.E(Double.NaN);
    }
    
    public IDataInfo processDataInfo(IDataInfo incomingDataInfo) {
        super.processDataInfo(incomingDataInfo);

        int n = incomingDataInfo.getLength();
        ratio = new DataDoubleArray(new int[]{n,n});
        ratioError = new DataDoubleArray(new int[]{n,n});

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

        String incomingLabel = incomingDataInfo.getLabel();
        DataDoubleArray.DataInfoDoubleArray ratioInfo = new DataDoubleArray.DataInfoDoubleArray(
                incomingLabel+" ratio", Null.DIMENSION, new int[]{n,n});
        DataDoubleArray.DataInfoDoubleArray ratioErrorInfo = new DataDoubleArray.DataInfoDoubleArray(
                incomingLabel+" ratio error", Null.DIMENSION, new int[]{n,n});

        ratioInfo.addTag(ratioTag);
        ratioErrorInfo.addTag(ratioErrorTag);
        
        subDataInfo = (IDataInfo[])Arrays.addObject(subDataInfo, ratioInfo);
        subDataInfo = (IDataInfo[])Arrays.addObject(subDataInfo, ratioErrorInfo);
        groupFactory.setSubDataInfo(subDataInfo);

        dataInfo = groupFactory.makeDataInfo();
        reset();
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
    protected DataDoubleArray ratio, ratioError;
    private final DataTag ratioTag, ratioErrorTag;
}
