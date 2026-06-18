/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSource;
import etomica.space.Vector;

/**
 * Meter for evaluation of the soft-potential pressure in a box.
 * Requires that temperature be set in order to calculation ideal-gas
 * contribution to pressure; default is to use zero temperature, which
 * causes this contribution to be omitted.
 *
 * @author David Kofke
 */

public class MeterReflected implements IDataSource {

    protected final CoordinateDefinition coordinateDefinition;
    protected final Box box, boxReflected;
    protected final IDataSource ds, dsr;

    public MeterReflected(IDataSource dataSource, IDataSource dataSourcer, Box box, Box boxReflected, CoordinateDefinition coordinateDefinition) {
        this.ds = dataSource;
        this.dsr = dataSourcer;
        this.box = box;
        this.boxReflected = boxReflected;
        this.coordinateDefinition = coordinateDefinition;
    }

    public IDataInfo getDataInfo() {
        return ds.getDataInfo();
    }

    public DataTag getTag() {
        return ds.getTag();
    }

    /**
     * Computes total pressure in box by summing virial over all pairs, and adding
     * ideal-gas contribution.
     */
    public IData getData() {
        int n = box.getLeafList().size();
        for (int i=0; i<n; i++) {
            IAtom a = box.getLeafList().get(i);
            IAtom ar = boxReflected.getLeafList().get(i);
            Vector dr = box.getSpace().makeVector();
            dr.Ev1Mv2(a.getPosition(), coordinateDefinition.getLatticePosition(a));
            ar.getPosition().Ev1Mv2(coordinateDefinition.getLatticePosition(a), dr);
        }
        IData d = ds.getData();
        double d2 = d.getValue(2);
        IData dr = dsr.getData();
        double d2r = dr.getValue(2);
        System.out.println(d2+" "+d2r+" "+(d2+d2r)/2);
        dr.PE(d);
        dr.TE(0.5);
        return dr;
    }

}
