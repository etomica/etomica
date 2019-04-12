/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations.KnottedPolymer;

import etomica.action.IAction;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.DataSourceUniform.LimitType;
import etomica.data.histogram.Histogram;
import etomica.data.histogram.HistogramSimple;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.math.DoubleRange;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Null;

/**
 * Meter for tabulation of the atomic radial distribution function (RDF).  The
 * meter takes data via actionPerformed and returns the average RDF via
 * getData.
 *
 * @author David Kofke
 */
public class MeterBondLength implements IAction, IDataSource, DataSourceIndependent, java.io.Serializable {

    private static final long serialVersionUID = 1L;
    protected final Space space;
    protected final DataSourceUniform xDataSource;
    protected final DataTag tag;
    protected final Vector dr;
    protected Box box;
    protected long[] gSum;
    protected DataFunction data;
    protected DataDoubleArray rData;
    protected double xMax;
    protected long callCount;
    protected AtomType type1, type2;
    private IDataInfo dataInfo;
    private String name;
    private final Histogram hist;
    private int f, l;

    /**
     * Creates meter with default to compute pair correlation for all
     * leaf atoms in a box.
     *
     * @param
     */
    public MeterBondLength(Box box, int f, int l) {
        this.space = box.getSpace();
        this.box = box;
        this.f = f;
        this.l = l;

        xDataSource = new DataSourceUniform("r", Length.DIMENSION);
        xDataSource.setTypeMax(LimitType.HALF_STEP);
        xDataSource.setTypeMin(LimitType.HALF_STEP);

        rData = (DataDoubleArray) xDataSource.getData();
        data = new DataFunction(new int[]{rData.getLength()});
        gSum = new long[rData.getLength()];
        dataInfo = new DataInfoFunction("g(r)", Null.DIMENSION, this);

        dr = space.makeVector();
        tag = new DataTag();
        dataInfo.addTag(tag);
        hist = new HistogramSimple(500, new DoubleRange(0, 10));
        xDataSource.setXMin(0);
        xDataSource.setXMax(hist.getXRange().maximum());
        xDataSource.setNValues(hist.getNBins());
        reset();
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    /**
     * Zero's out the RDF sum tracked by this meter.
     */
    public void reset() {
        data = new DataFunction(new int[]{hist.getNBins()}, hist.getHistogram());
        dataInfo = new DataInfoFunction("Bond length", Length.DIMENSION, this);
        dataInfo.addTag(tag);
        hist.reset();
    }

    /**
     * Takes the RDF for the current configuration of the given box.
     */
    public void actionPerformed() {

        IAtomList atomList = box.getLeafList();
        for (int i = 0; i < f; i++) {
            IAtom prev = atomList.get(0);
            for (int j = 0; j < l; j++) {
                int k = 1 + i * l + j;
                IAtom a = atomList.get(k);
                double r2 = prev.getPosition().Mv1Squared(a.getPosition());
//                if (r2 > 1.5*1.5) {
//                    throw new RuntimeException("oops "+i+" "+j+" "+k+" "+Math.sqrt(r2)+" "+prev.getPosition()+" "+a.getPosition());
//                }
                prev = a;
//                if (j==0) {
////                    continue;
//                }else{
//                    continue;
//                }
                hist.addValue(Math.sqrt(r2));


            }
        }
    }

    /**
     * Returns the RDF, averaged over the calls to actionPerformed since the
     * meter was reset or had some parameter changed (xMax or # of bins).
     */
    public IData getData() {
        hist.getHistogram();
        return data;
    }

    public DataSourceUniform getXDataSource() {
        return xDataSource;
    }

    public DataDoubleArray getIndependentData(int i) {
        return (DataDoubleArray) xDataSource.getData();
    }

    public DataInfoDoubleArray getIndependentDataInfo(int i) {
        return (DataInfoDoubleArray) xDataSource.getDataInfo();
    }

    public DataTag getIndependentTag() {
        return xDataSource.getTag();
    }

    public int getIndependentArrayDimension() {
        return 1;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }
}
