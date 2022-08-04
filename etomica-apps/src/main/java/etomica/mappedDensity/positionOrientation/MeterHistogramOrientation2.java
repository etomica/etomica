/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.mappedDensity.positionOrientation;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.molecule.IMolecule;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

public class MeterHistogramOrientation2 implements IDataSource, DataSourceIndependent {

    protected DataFunction data;
    protected DataFunction.DataInfoFunction dataInfo;
    protected DataTag tag;
    protected DataSourceUniform cosSource, perpSource;
    protected double[] perps, widths;
    protected final Box box;
    protected int d;
    protected double S;

    public MeterHistogramOrientation2(Box box, int d, int nBinsCos, double[] perps, double[] widths) {
        this.box = box;
        this.d = d;
        int D = box.getSpace().D();
        this.perps = perps;
        this.widths = widths;
        if(D==2) {
            cosSource = new DataSourceUniform("angle", Null.DIMENSION, nBinsCos, -Math.PI, Math.PI, DataSourceUniform.LimitType.HALF_STEP, DataSourceUniform.LimitType.HALF_STEP);
            S = box.getBoundary().getBoxSize().getX(0);
        }
        else {
            cosSource = new DataSourceUniform("cosine", Null.DIMENSION, nBinsCos, 0, 1, DataSourceUniform.LimitType.HALF_STEP, DataSourceUniform.LimitType.HALF_STEP);
            S = box.getBoundary().getBoxSize().getX(0) * box.getBoundary().getBoxSize().getX(1);
        }
        // not useful, but need to provide something as independent data
        perpSource = new DataSourceUniform("perp", Null.DIMENSION, perps.length*widths.length, 0, 1);
        data = new DataFunction(new int[]{widths.length*perps.length,nBinsCos});
        dataInfo = new DataFunction.DataInfoFunction("histogram", Null.DIMENSION, this);
        tag = new DataTag();
        dataInfo.addTag(tag);
    }

    @Override
    public IData getData() {
        data.E(0);
        double[] y = data.getData();
        int D = box.getSpace().D();
        for (IMolecule molecule : box.getMoleculeList()) {
            Vector dr = box.getSpace().makeVector();
            IAtomList atoms = molecule.getChildList();
            dr.Ev1Mv2(atoms.get(1).getPosition(), atoms.get(0).getPosition());
            box.getBoundary().nearestImage(dr);
            Vector xyz = box.getSpace().makeVector();
            xyz.E(atoms.get(0).getPosition());
            xyz.PEa1Tv1(0.5, dr);//location of center of dimer

            double dz = Math.abs(xyz.getX(d));
            for (int iw=0; iw<widths.length; iw++) {
                double inc = 1/(widths[iw] * S) * (cosSource.getNValues() / (2. * Math.PI));
                for (int ip=0; ip<perps.length; ip++) {
                    double iz = perps[ip];
                    if (Math.abs(dz - iz) < widths[iw] / 2) {
                        double drd;
                        if (D == 2) drd = Math.atan2(dr.getX(1), dr.getX(0));// /(Math.PI);
                        else throw new RuntimeException("not set up for 3D in MeterHistogramOrientation");
                        int cosIdx = cosSource.getIndex(drd);
                        y[(iw * perps.length + ip) * cosSource.getNValues() + cosIdx] += inc;
                    }
                }
            }
        }
        return data;
    }

    @Override
    public DataTag getTag() {
        return tag;
    }

    @Override
    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    @Override
    public DataDoubleArray getIndependentData(int i) {
        return (DataDoubleArray) (i==0 ? perpSource.getData() : cosSource.getData());
    }

    @Override
    public DataDoubleArray.DataInfoDoubleArray getIndependentDataInfo(int i) {
        return (DataDoubleArray.DataInfoDoubleArray) (i==0 ? perpSource.getDataInfo() : cosSource.getDataInfo());
    }

    @Override
    public int getIndependentArrayDimension() {
        return 2;
    }

    @Override
    public DataTag getIndependentTag() {
        return null;
    }
}
