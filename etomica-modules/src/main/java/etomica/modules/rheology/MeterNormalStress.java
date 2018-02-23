/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.rheology;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.data.DataSourceScalar;
import etomica.space.Space;
import etomica.units.dimensions.Null;

/**
 * Meter to measure normal stress of a polymer in a shear field.
 *
 * @author Andrew Schultz
 */
public class MeterNormalStress extends DataSourceScalar {

    public MeterNormalStress(Space space) {
        super("normal stress coefficient", Null.DIMENSION);
        dr = space.makeVector();
    }

    public void setIntegrator(IntegratorPolymer newIntegrator) {
        integrator = newIntegrator;
    }

    public void setBox(Box newBox) {
        box = newBox;
    }

    public double getDataAsScalar() {
        double shearRate = integrator.getShearRate();
        if (shearRate == 0) {
            return Double.NaN;
        }
        double b = integrator.getB();
        IAtomList list = box.getMoleculeList().get(0).getChildList();
        double s = 0;
        for (int i = 0; i<list.size()-1; i++) {
            Vector p0 = list.get(i).getPosition();
            Vector p1 = list.get(i+1).getPosition();
            dr.Ev1Mv2(p1, p0);
            double fQ = 1+b*dr.squared();
            double dr0 = dr.getX(d[0])*fQ;
            double dr1 = dr.getX(d[1])*fQ;
            s += (dr0*dr0 - dr1*dr1);
        }
        s /= shearRate;
        if (doDouble) {
            s /= shearRate;
        }
        return s;
    }
    
    public void setDoDouble(boolean newDoDouble) {
        doDouble = newDoDouble;
    }

    public void setDims(int[] newDims) {
        d = newDims;
    }

    private static final long serialVersionUID = 1L;
    protected Box box;
    protected Vector dr;
    protected IntegratorPolymer integrator;
    protected int[] d;
    protected boolean doDouble;
}
