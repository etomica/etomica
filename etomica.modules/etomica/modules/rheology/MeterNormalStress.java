/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.rheology;

import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.data.DataSourceScalar;
import etomica.space.ISpace;
import etomica.units.Null;

/**
 * Meter to measure normal stress of a polymer in a shear field.
 *
 * @author Andrew Schultz
 */
public class MeterNormalStress extends DataSourceScalar {

    public MeterNormalStress(ISpace space) {
        super("normal stress coefficient", Null.DIMENSION);
        dr = space.makeVector();
    }

    public void setIntegrator(IntegratorPolymer newIntegrator) {
        integrator = newIntegrator;
    }

    public void setBox(IBox newBox) {
        box = newBox;
    }

    public double getDataAsScalar() {
        double shearRate = integrator.getShearRate();
        if (shearRate == 0) {
            return Double.NaN;
        }
        double b = integrator.getB();
        IAtomList list = box.getMoleculeList().getMolecule(0).getChildList();
        double s = 0;
        for (int i=0; i<list.getAtomCount()-1; i++) {
            IVector p0 = list.getAtom(i).getPosition();
            IVector p1 = list.getAtom(i+1).getPosition();
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
    protected IBox box;
    protected IVectorMutable dr;
    protected IntegratorPolymer integrator;
    protected int[] d;
    protected boolean doDouble;
}
