/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.integrator.mcmove.MCMoveOverlapListenerFasterer;
import etomica.units.dimensions.Null;
import etomica.units.dimensions.Quantity;

/**
 * Free energy data source that computes free energy differences based on
 * overlap sampling performed by an MCMoveOverlapListener.
 * 
 * @author Andrew Schultz
 */
public class DataSourceFEFasterer implements IDataSource, DataSourceIndependent {

    protected final MCMoveOverlapListenerFasterer mcMoveOverlapMeter;
    protected DataFunction data;
    protected DataInfoFunction dataInfo;
    protected DataDoubleArray xData;
    protected DataInfoDoubleArray xDataInfo;
    protected final DataTag tag, xTag;
    protected boolean doSubtractComb, doAvgDef;

    public DataSourceFEFasterer(MCMoveOverlapListenerFasterer mcMoveOverlapMeter) {
        tag = new DataTag();
        xTag = new DataTag();
        this.mcMoveOverlapMeter = mcMoveOverlapMeter;
        xData = new DataDoubleArray(0);
        xDataInfo = new DataInfoDoubleArray("bar", Null.DIMENSION, new int[]{0});
        xDataInfo.addTag(xTag);
        data = new DataFunction(new int[]{0});
        dataInfo = new DataInfoFunction("foo", Null.DIMENSION, this);
        dataInfo.addTag(tag);
    }
    
    /**
     * If true, the free energy differences will not include the combinatoric
     * contributions (they will be the free energy of defect formation).
     * 
     * http://dx.doi.org/10.1073/pnas.1211784109
     */
    public void setSubtractComb(boolean doSubtractComb) {
        this.doSubtractComb = doSubtractComb;
    }

    public void setAvgDef(boolean doAvgDef) {
        this.doAvgDef = doAvgDef;
    }

    public IData getData() {
        double[] ratios = mcMoveOverlapMeter.getRatios().clone();
        if (ratios == null || ratios.length==0) return data;
        if (ratios.length != dataInfo.getLength()-1) {
            getDataInfo();
        }
        double[] y = data.getData();
        long[] numInsert = mcMoveOverlapMeter.getNumInsert();
        long[] numDelete = mcMoveOverlapMeter.getNumDelete();
        int n0 = mcMoveOverlapMeter.getMinNumAtoms();
        double wTot = 0, daDefTot = 0;
        for (int i=0; i<ratios.length; i++) {
            if (!Double.isNaN(ratios[i])) {
                y[i] = Math.log(ratios[i]);
                if (doSubtractComb) {
                    y[i] += Math.log((n0+i+1)/(ratios.length-i));
                    if (doAvgDef) {
                        double w = 1.0/(1.0/numInsert[i]+1.0/numDelete[i+1]);
                        daDefTot += y[i]*w;
                        wTot += w;
                    }
                }
            }
        }
        if (doAvgDef) {
            double daDef = daDefTot/wTot;
            for (int i=0; i<ratios.length; i++) {
                y[i] = daDef;
            }
        }
        return data;
    }

    public DataTag getTag() {
        return tag;
    }

    public IDataInfo getDataInfo() {
        double[] ratios = mcMoveOverlapMeter.getRatios();
        if (ratios == null) return dataInfo;
        if (ratios.length != dataInfo.getLength()) {
            xDataInfo = new DataInfoDoubleArray("N", Quantity.DIMENSION, new int[]{ratios.length});
            xDataInfo.addTag(tag);
            xData = new DataDoubleArray(ratios.length);
            double[] x = xData.getData();
            int n0 = mcMoveOverlapMeter.getMinNumAtoms();
            for (int i=0; i<ratios.length; i++) {
                x[i] = n0+i;
            }
            
            dataInfo = new DataInfoFunction("FE", Null.DIMENSION, this);
            dataInfo.addTag(tag);
            data = new DataFunction(new int[]{ratios.length});
        }
        return dataInfo;
    }

    public DataDoubleArray getIndependentData(int i) {
        return xData;
    }

    public DataInfoDoubleArray getIndependentDataInfo(int i) {
        return xDataInfo;
    }

    public int getIndependentArrayDimension() {
        return 1;
    }

    public DataTag getIndependentTag() {
        return xTag;
    }

}
