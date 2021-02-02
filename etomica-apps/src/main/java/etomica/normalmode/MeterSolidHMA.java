/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.potential.compute.PotentialCallback;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

public class MeterSolidHMA implements IDataSource, PotentialCallback {

    protected final int dim;
    protected final DataTag tag;
    protected DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected final PotentialCompute potentialMaster;
    protected double temperature;
    protected double latticeEnergy, latticePressure;
    protected final Box box;
    protected double pRes, bpHarm;
    protected final boolean doD2;
    protected final CoordinateDefinition coordinteDefinition;
    protected final Vector dr;
    protected double d2sum;

    public MeterSolidHMA(Space space, PotentialCompute potentialCompute, CoordinateDefinition coordinateDefinition, boolean doD2) {
        this.coordinteDefinition = coordinateDefinition;
        tag = new DataTag();
        this.potentialMaster = potentialCompute;
        dim = space.D();
        box = coordinateDefinition.getBox();
        latticeEnergy = potentialCompute.computeAll(false);

        latticePressure = -potentialCompute.getLastVirial() / (box.getBoundary().volume() * dim);

        int n = doD2 ? 7 : 5;
        dataInfo = new DataInfoDoubleArray("Stuff", Null.DIMENSION, new int[]{n});
        dataInfo.addTag(tag);
        data = new DataDoubleArray(n);
        this.doD2 = doD2;
        dr = space.makeVector();
    }

    @Override
    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    @Override
    public DataTag getTag() {
        return tag;
    }

    public void setTemperature(double temperature) {
        this.temperature = temperature;
        bpHarm = pRes / temperature;
    }

    public void setPRes(double pRes) {
        this.pRes = pRes;
        bpHarm = pRes / temperature;
    }

    /**
     * Computes total pressure in box by summing virial over all pairs, and adding
     * ideal-gas contribution.
     */
    @Override
    public IData getData() {
        PotentialCallback callback = doD2 ? this : null;
        d2sum = 0;
        double uSum = potentialMaster.computeAll(true, callback);
        double[] x = data.getData();
        double V = box.getBoundary().volume();
        double rho = box.getMoleculeList().size() / V;
        double virial = potentialMaster.getLastVirial();
        double measuredP = temperature * rho - virial / (dim * V);
        int N = box.getMoleculeList().size();
        x[0] = uSum / N;
        x[1] = measuredP;

        Vector[] forces = potentialMaster.getForces();
        double fdr = 0;
        IAtomList atoms = box.getLeafList();

        for (IAtom a : atoms) {
            Vector f = forces[a.getLeafIndex()];
            dr.Ev1Mv2(a.getPosition(), coordinteDefinition.getLatticePosition(a));
            box.getBoundary().nearestImage(dr);
            fdr += f.dot(dr);
        }

        double buc = (0.5 * fdr + (uSum - latticeEnergy)) / temperature / N;
        x[2] = buc;
        double vol = box.getBoundary().volume();
        // P = Plat + Pres + x[5]
        double density = N / vol;
        double fV = (bpHarm - N / V) / (dim * (N - 1));
        double Zc = (-virial / (dim * V) + fV * fdr - latticePressure) / (density * temperature);

        x[3] = Zc;
        // Pc = x[5]
        // Zc = x[5] / (rho*T)
        // this is dbAc/dv2 at constant Y (for LJ)


        x[4] = (4 * buc - Zc) * density * density / 2;
        //x[4] = (N-1)*1.5/N + (0.5 * pc.getDADBSum() + uSum) / temperature / N;

        if (doD2) {
            double y = latticeEnergy + dr.getD() * (N - 1) * temperature / 2;
            x[5] = (uSum - y) * (uSum - y);
//            System.out.println(uSum+" "+y+" "+x[5]);

            x[6] = -0.25 * (fdr + d2sum) / temperature + (x[2] * x[2]) * N * N;
            //x[7] = -0.25 * (pc.getDADBSum() + pc.getD2Sum()) / temperature + x[4]*x[4]*N*N;
//            x[7] = pc.getD2Sum();
        }

        return data;
    }

    @Override
    public void pairCompute(int iAtom, int jAtom, Vector dr, double[] u012) {
        double d2u = u012[2];
        double r2 = dr.squared();
        double dud2u = (u012[1] - d2u) / (r2 * r2);
        int D = dr.getD();
        Vector dri = Vector.d(D);
        Vector drj = Vector.d(D);
        IAtom ia = box.getLeafList().get(iAtom);
        IAtom ja = box.getLeafList().get(jAtom);
        dri.Ev1Mv2(ia.getPosition(), coordinteDefinition.getLatticePosition(ia));
        box.getBoundary().nearestImage(dri);
        drj.Ev1Mv2(ja.getPosition(), coordinteDefinition.getLatticePosition(ja));
        box.getBoundary().nearestImage(drj);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                double der2 = dr.getX(i) * dr.getX(j) * dud2u;
                if (i == j) der2 -= u012[1] / r2;
                d2sum += der2 * (2 * dri.getX(i) * drj.getX(j) - dri.getX(i) * dri.getX(j) - drj.getX(i) * drj.getX(j));
            }
        }
    }
}
