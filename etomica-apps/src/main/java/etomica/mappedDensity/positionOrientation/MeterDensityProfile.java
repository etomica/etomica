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

/**
 * Measures histogram of molecule positions in z direction, regardless of orientation
 * Does not exploit symmetry about centerline between walls
 */
public class MeterDensityProfile implements IDataSource, DataSourceIndependent {

    protected DataFunction data;
    protected DataFunction.DataInfoFunction dataInfo;
    protected DataTag tag;
    protected DataSourceUniform perpSource;
    protected final Box box;
    protected final int d;
    protected double S;

    public MeterDensityProfile(Box box, int d, int nBinsPerp) {
        this.box = box;
        this.d = d;
        double L = box.getBoundary().getBoxSize().getX(d);
        perpSource = new DataSourceUniform("perp", Length.DIMENSION, nBinsPerp, 0, L, DataSourceUniform.LimitType.HALF_STEP, DataSourceUniform.LimitType.HALF_STEP);
        data = new DataFunction(new int[]{nBinsPerp});
        dataInfo = new DataFunction.DataInfoFunction("histogram", Null.DIMENSION, this);
        tag = new DataTag();
        dataInfo.addTag(tag);
        S = box.getBoundary().getBoxSize().getX(0);  //area
        if(d == 3) S *= box.getBoundary().getBoxSize().getX(1);
    }

    @Override
    public IData getData() {
        data.E(0);
        double[] y = data.getData();
        double inc = perpSource.getNValues() / (S * box.getBoundary().getBoxSize().getX(d));
        double Lz = box.getBoundary().getBoxSize().getX(d);
        for (IMolecule molecule : box.getMoleculeList()) {
            Vector dr = box.getSpace().makeVector();
            IAtomList atoms = molecule.getChildList();
            dr.Ev1Mv2(atoms.get(1).getPosition(), atoms.get(0).getPosition());
            box.getBoundary().nearestImage(dr);
            Vector xyz = box.getSpace().makeVector();
            xyz.E(atoms.get(0).getPosition());
            xyz.PEa1Tv1(0.5, dr);
            double dz = xyz.getX(d) + Lz/2;
            int perpIdx = perpSource.getIndex(dz);
            y[perpIdx] += inc;
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
        return (DataDoubleArray) perpSource.getData();
    }

    @Override
    public DataDoubleArray.DataInfoDoubleArray getIndependentDataInfo(int i) {
        return (DataDoubleArray.DataInfoDoubleArray) perpSource.getDataInfo();
    }

    @Override
    public int getIndependentArrayDimension() {
        return 1;
    }

    @Override
    public DataTag getIndependentTag() {
        return null;
    }
}
