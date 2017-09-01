/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.data.DataSourceIndependent;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.integrator.mcmove.MCMoveInsertDeleteBiased;
import etomica.integrator.mcmove.MCMoveOverlapListener;
import etomica.units.dimensions.Null;
import etomica.units.dimensions.Quantity;

/**
 * Like DataSourceFEHistogram, returns a probability histogram for N based on
 * free energies.  However, this class does a weighted average to compute the
 * average defect free energy (weighted by how much each N is visited) and then
 * assumes that the defect free energy is constant for all N.  This yields a
 * fairly accurate description fairly quickly and tends to be more accurate for
 * undersampled N values.
 * 
 * @author Andrew Schultz
 */
public class DataSourceImposedFEFit implements IEtomicaDataSource, DataSourceIndependent {

    protected final MCMoveOverlapListener mcMoveOverlapMeter;
    protected DataFunction data;
    protected DataInfoFunction dataInfo;
    protected DataDoubleArray xData;
    protected DataInfoDoubleArray xDataInfo;
    protected final DataTag tag, xTag;
    protected double bmu;
    protected final MCMoveInsertDeleteBiased mcMoveID;
    
    public DataSourceImposedFEFit(MCMoveOverlapListener mcMoveOverlapMeter, MCMoveInsertDeleteBiased mcMoveID, double bmu) {
        tag = new DataTag();
        xTag = new DataTag();
        this.mcMoveOverlapMeter = mcMoveOverlapMeter;
        this.mcMoveID = mcMoveID;
        xDataInfo = new DataInfoDoubleArray("bar", Null.DIMENSION, new int[]{0});
        data = new DataFunction(new int[]{0});
        dataInfo = new DataInfoFunction("foo", Null.DIMENSION, this);
        this.bmu = bmu;
    }
    
    public void setMu(double newBMu) {
        bmu = newBMu;
    }

    public IData getData() {
        double[] ratios = mcMoveOverlapMeter.getRatios().clone();
        if (ratios == null || ratios.length==0) return data;
        if (ratios.length != dataInfo.getLength()-1) {
            getDataInfo();
        }
        long[] numInsert = mcMoveOverlapMeter.getNumInsert();
        long[] numDelete = mcMoveOverlapMeter.getNumDelete();
        int n0 = mcMoveOverlapMeter.getMinNumAtoms();
        double wTot = 0, daDefTot = 0;
        for (int i=0; i<ratios.length; i++) {
            if (!Double.isNaN(ratios[i])) {
                double iDaDef = Math.log(ratios[i]);
                iDaDef += Math.log((n0+i+1)/(ratios.length-i));
                double w = 1.0/(1.0/numInsert[i]+1.0/numDelete[i+1]);
                daDefTot += iDaDef*w;
                wTot += w;
            }
        }
        double daDef = daDefTot/wTot;

        double p = 1;
        double tot = 0;
        for (int i=ratios.length-1; i>=0; i--) {
            tot += p;
            p /= Math.exp(daDef+bmu)/(n0+i+1)*(ratios.length-i);
        }
        tot += p;
        double[] y = data.getData();
        double p2 = 1;
        for (int i=ratios.length; i>=0; i--) {
            y[i] = p2 == 0 ? Double.NaN : p2/tot;
            if (i==0) break;
            p2 /= Math.exp(daDef+bmu)/(n0+i)*(ratios.length-i+1);
        }
        return data;
    }

    public DataTag getTag() {
        return tag;
    }

    public IEtomicaDataInfo getDataInfo() {
        double[] ratios = mcMoveOverlapMeter.getRatios();
        if (ratios == null) return dataInfo;
        if (ratios.length != dataInfo.getLength()-1) {
            xDataInfo = new DataInfoDoubleArray("N", Quantity.DIMENSION, new int[]{ratios.length+1});
            xDataInfo.addTag(tag);
            xData = new DataDoubleArray(ratios.length+1);
            double[] x = xData.getData();
            int n0 = mcMoveOverlapMeter.getMinNumAtoms();
            for (int i=0; i<=ratios.length; i++) {
                x[i] = n0+i;
            }
            
            dataInfo = new DataInfoFunction("FE Histogram", Null.DIMENSION, this);
            dataInfo.addTag(tag);
            data = new DataFunction(new int[]{ratios.length+1});
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
