/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.mappedvirial;

import etomica.action.IAction;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.histogram.Histogram;
import etomica.data.histogram.HistogramNotSoSimple;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.math.DoubleRange;
import etomica.potential.IPotential2;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Vector;
import etomica.units.dimensions.Force;
import etomica.units.dimensions.Length;

public class MeterMeanForce implements IDataSource, DataSourceIndependent, IAction {

    protected final PotentialCompute potentialMaster;
    protected final Box box;
    protected final IPotential2 p2;
    protected final Vector dr, fij;
    protected final DataFunction data;
    protected final DataInfoFunction dataInfo;
    protected final DataTag tag, xTag;
    protected final HistogramNotSoSimple hist, hist2;
    protected final DataDoubleArray xData;
    protected final DataInfoDoubleArray xDataInfo;

    public MeterMeanForce(PotentialCompute potentialMaster, IPotential2 p2, Box box, int nbins) {
        this.p2 = p2;
        this.box = box;
        this.potentialMaster = potentialMaster;
        dr = box.getSpace().makeVector();
        fij = box.getSpace().makeVector();

        xData = new DataDoubleArray(nbins);
        xDataInfo = new DataInfoDoubleArray("r", Length.DIMENSION, new int[]{nbins});
        xTag = new DataTag();
        xDataInfo.addTag(xTag);
        
        data = new DataFunction(new int[]{nbins});
        dataInfo = new DataInfoFunction("mean force", Force.DIMENSION, this);
        tag = new DataTag();
        dataInfo.addTag(tag);
        
        hist = new HistogramNotSoSimple(nbins, new DoubleRange(0,p2.getRange()));
        hist2 = new HistogramNotSoSimple(nbins, new DoubleRange(0,p2.getRange()));
    }
    
    public Histogram getHistogram() {
        return hist;
    }
    
    public Histogram getHistogram2() {
        return hist2;
    }

    public void actionPerformed() {
    }
    
    public IData getData() {
        potentialMaster.computeAll(true);
        Vector[] forces = potentialMaster.getForces();
        IAtomList list = box.getLeafList();

        int n = list.size();
        double rc2 = p2.getRange()*p2.getRange();
        for (int i=0; i<n; i++) {
            IAtom a = list.get(i);
            Vector fi = forces[i];
            for (int j=i+1; j<n; j++) {
                IAtom b = list.get(j);
                dr.Ev1Mv2(b.getPosition(),a.getPosition());
                box.getBoundary().nearestImage(dr);
                double r2 = dr.squared();
                if (r2 > rc2) continue;
                double r = Math.sqrt(r2);
                Vector fj = forces[j];
                fij.Ev1Mv2(fj,fi);
                double fdr = 0.5*fij.dot(dr)/r;
                hist.addValue(r, fdr);
                hist2.addValue(r, fdr*fdr);
            }
        }


        double[] h = hist.getHistogram();
        double[] y = data.getData();
        System.arraycopy(h, 0, y, 0, h.length);
        return data;
    }

    public DataTag getTag() {
        return tag;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataDoubleArray getIndependentData(int i) {
        double[] r = hist.xValues();
        double[] x = xData.getData();
        System.arraycopy(r, 0, x, 0, r.length);
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
