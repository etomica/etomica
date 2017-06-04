/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.droplet;

import etomica.atom.IAtomList;
import etomica.space.Vector;
import etomica.box.Box;
import etomica.atom.AtomFilter;
import etomica.atom.AtomFilterStatic;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.units.Null;

public class MeterDeformation implements IEtomicaDataSource {

    public MeterDeformation(Space space) {
        data = new DataDoubleArray(2);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("deformation", Null.DIMENSION, new int[]{2});
        tag = new DataTag();
        f = new double[10];
        for (int i=0; i<f.length; i++) {
            f[i] = 1 + (i+1)/3.0;
        }

        center = space.makeVector();
        dr = space.makeVector();
        rg = space.makeVector();
        moment = space.makeTensor();
        workTensor = space.makeTensor();
        setFilter(AtomFilterStatic.ACCEPT_ALL);
    }

    public void setBox(Box newBox) {
        box = newBox;
    }

    public Box getBox() {
        return box;
    }

    public void setFilter(AtomFilter newFilter) {
        filter = newFilter;
    }

    public AtomFilter getFilter() {
        return filter;
    }

    public IData getData() {
        IAtomList leafList = box.getLeafList();
        center.E(0);
        for (int i=0; i<leafList.getAtomCount(); i++) {
            if (!filter.accept(leafList.getAtom(i))) {
                continue;
            }
            center.PE(leafList.getAtom(i).getPosition());
        }
        center.TE(1.0/leafList.getAtomCount());

        rg.E(0);
        moment.E(0);
        for (int i=0; i<leafList.getAtomCount(); i++) {
            if (!filter.accept(leafList.getAtom(i))) {
                continue;
            }
            dr.Ev1Mv2(leafList.getAtom(i).getPosition(), center);
            workTensor.Ev1v2(dr, dr);
            moment.PE(workTensor);

        }
        
        moment.TE(1.0/leafList.getAtomCount());
        rg.setX(0, moment.component(0,0));
        rg.setX(1, moment.component(1,1));
        rg.setX(2, moment.component(2,2));
        // ugh
        if (rg.getX(0) > rg.getX(1)) {
            double t = rg.getX(0);
            rg.setX(0, rg.getX(1));
            rg.setX(1, t);
        }
        if (rg.getX(2) < rg.getX(0)) {
            double t = rg.getX(0);
            rg.setX(0, rg.getX(2));
            rg.setX(2, rg.getX(1));
            rg.setX(1, t);
        }
        else if (rg.getX(2) < rg.getX(1)) {
            double t = rg.getX(1);
            rg.setX(1, rg.getX(2));
            rg.setX(2, t);
        }
        // now rg0 <= rg1 <= rg2

        double b1 = moment.trace();

        double b2 = 0;
        for (int i=0; i<moment.D(); i++) {
            for (int j=0; j<moment.D(); j++) {
                double m = moment.component(i,j);
                b2 += m*m;
            }
        }

        double b3 = moment.determinant();

        double c1 = b1 / (3.0*Math.pow(b3, 1.0/3.0)) - 1.0;
//        double c2 = b2 / (3.0*Math.pow(b3, 2.0/3.0)) - 1.0;

        // rg1 is closer to rg0 or rg2?
//        double delta = c2 - 4.0*c1;
//        double factor = Math.signum(delta);
        double factor = 1;  // assume prolate
        if (rg.getX(1)*rg.getX(1) > rg.getX(0)*rg.getX(2)) {
            factor = -1;    // oblate
        }

        double eOld = 0;
        double eOldOld = 0;
        double e = factor * Math.sqrt(c1/3.0);
        int iter = 0;
        while (iter < 20 && Math.abs(e-eOld) > 1.e-9 && (Math.abs(e-eOld) <= Math.abs(e-eOldOld) || (e-eOld)*(e-eOldOld) < 0)) {
            eOldOld = eOld;
            eOld = e;
            double denominator = 1.0;
            for (int i=0; i<10; i++) {
                denominator += f[i] * Math.pow(eOld, (i+1));
            }
            e = factor * Math.sqrt(c1 / denominator);
            iter++;
        }

        double[] x = data.getData();
        x[0] = Math.pow(125*b3, 1.0/6.0);
        if (e > 0.02) {
            double a = Math.pow(1-e, 1.5);
            x[1] = (1.0 - a) / (1 + a);
        }
        else {
            double e2 = e*e;
            x[1] = 0.75*e + 0.375*e2 + 7.0/64.0*e2*e - 3.0/128.0*e2*e2 + 33.0/512.0*e2*e2*e - 61.0/1024.0*e2*e2*e2; 
        }
        return data;
    }

    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    protected Box box;
    protected final DataDoubleArray data;
    protected final DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected final double[] f;
    protected final Vector center;
    protected final Vector dr;
    protected final Vector rg;
    protected final Tensor moment, workTensor;
    protected AtomFilter filter;
}
