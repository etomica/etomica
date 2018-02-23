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
 * Meter to measure the viscosity of a polymer in a shear field.
 *
 * @author Andrew Schultz
 */
public class MeterViscosity extends DataSourceScalar {

    public MeterViscosity(Space space) {
        super("viscosity", Null.DIMENSION);
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
        IAtomList list = box.getMoleculeList().getMolecule(0).getChildList();
        double v = 0;
        for (int i = 0; i<list.size()-1; i++) {
            Vector p0 = list.get(i).getPosition();
            Vector p1 = list.get(i+1).getPosition();
            dr.Ev1Mv2(p1, p0);
            v += dr.getX(0)*dr.getX(1)/(1+b*dr.squared());
        }
        return v/shearRate;
    }

    private static final long serialVersionUID = 1L;
    protected Box box;
    protected Vector dr;
    protected IntegratorPolymer integrator;
}
