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
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Null;

public class MeterHistogramOrientation implements IDataSource, DataSourceIndependent {

    protected DataFunction data;
    protected DataFunction.DataInfoFunction dataInfo;
    protected DataTag tag;
    protected DataSourceUniform perpSource, cosSource;
    protected final Box box;
    protected int d;

    public MeterHistogramOrientation(Box box, int d, int nBinsPerp, int nBinsCos) {
        this.box = box;
        this.d = d;
        double L = box.getBoundary().getBoxSize().getX(d);
        int D = box.getSpace().D();
        perpSource = new DataSourceUniform("perp", Length.DIMENSION, nBinsPerp, 0, L/2, DataSourceUniform.LimitType.HALF_STEP, DataSourceUniform.LimitType.HALF_STEP);
        cosSource = new DataSourceUniform(D==2 ? "angle" : "cosine", Null.DIMENSION, nBinsCos, 0, D==2 ? Math.PI/2 : 1, DataSourceUniform.LimitType.HALF_STEP, DataSourceUniform.LimitType.HALF_STEP);
        data = new DataFunction(new int[]{nBinsPerp,nBinsCos});
        dataInfo = new DataFunction.DataInfoFunction("histogram", Null.DIMENSION, this);
        tag = new DataTag();
        dataInfo.addTag(tag);
    }

    @Override
    public IData getData() {
        data.E(0);
        double[] y = data.getData();
        double inc = (perpSource.getNValues() / box.getBoundary().getBoxSize().getX(d)) * cosSource.getNValues();
        double Lz = box.getBoundary().getBoxSize().getX(d);
        int D = box.getSpace().D();
        for (IMolecule molecule : box.getMoleculeList()) {
            Vector dr = box.getSpace().makeVector();
            IAtomList atoms = molecule.getChildList();
            dr.Ev1Mv2(atoms.get(1).getPosition(), atoms.get(0).getPosition());
            box.getBoundary().nearestImage(dr);
            Vector xyz = box.getSpace().makeVector();
            xyz.E(atoms.get(0).getPosition());
            xyz.PEa1Tv1(0.5, dr);
            double dz = (0.5*Lz - Math.abs(xyz.getX(d)));
            int perpIdx = perpSource.getIndex(dz);
            dr.normalize();
            double drd = Math.abs(dr.getX(d));
            if (D == 2) drd = Math.acos(drd);
            int cosIdx = cosSource.getIndex(drd);
            y[perpIdx* cosSource.getNValues() + cosIdx] += inc;
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
