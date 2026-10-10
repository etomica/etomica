/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.potential.PotentialMaster;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

public class MeterSolidHMAFD implements IDataSource {

    protected Box box;
    protected Vector[] latticePositions;
    protected final DataDoubleArray data;
    protected final DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected final Space space;
    protected final CoordinateDefinition coordinateDefinition;
    protected final DataSourceScalar meterPE;
    protected final Vector dr;
    protected double V, rho;
    protected int N, nCells;
    protected double gV, gVV;
    protected double mV;
    protected Vector dri, drj, drij;
    protected final double temperature;
    protected final PotentialCompute potentialMaster;
    protected Vector[] rOrig;
    protected int n;

    public MeterSolidHMAFD(Space space, DataSourceScalar meterPE, PotentialMaster potentialMaster, CoordinateDefinition coordinateDefinition, double temperature, double[] elasticParams) {
        int nData = 5;
        data = new DataDoubleArray(nData);
        dataInfo = new DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{nData});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.space = space;
        this.box = coordinateDefinition.getBox();
        this.coordinateDefinition = coordinateDefinition;
        this.meterPE = meterPE;
        this.potentialMaster = potentialMaster;
        dr = space.makeVector();
        this.temperature = temperature;
        V = box.getBoundary().volume();
        N = box.getLeafList().size();
        rho = N / V;

        gV   = elasticParams[0];
        gVV  = elasticParams[1];
        mV = gV - 1.0 / 3.0;

        dri = space.makeVector();
        drj = space.makeVector();
        drij = space.makeVector();
        rOrig = box.getSpace().makeVectorArray(box.getLeafList().size());

        latticePositions = box.getSpace().makeVectorArray(N);
        for (int i = 0; i < N; i++) {
            latticePositions[i].E(coordinateDefinition.getLatticePosition(box.getLeafList().get(i)));
        }

        n = 0;
    }


    public IData getData() {
        n++;
        double[] x = data.getData();

        for (IAtom atom : box.getLeafList()) {
            rOrig[atom.getLeafIndex()].E(atom.getPosition());
        }

        double Fdr = 0;
        Box box = coordinateDefinition.getBox();
        potentialMaster.computeAll(true);

        double uSum = potentialMaster.getLastEnergy();
        double virial = potentialMaster.getLastVirial();
        Vector[] forces = potentialMaster.getForces();

        for (IAtom atom : box.getLeafList()) {
            int iAtom = atom.getLeafIndex();
            Vector F = forces[iAtom];
            dri.Ev1Mv2(atom.getPosition(), latticePositions[iAtom]);
            box.getBoundary().nearestImage(dri);
            Fdr += F.dot(dri);
        }

        //Uc
        x[0] = uSum / N;
        //Um
        x[1] = uSum / N + 3.0 * (N - 1.0) / 2.0 / N * temperature  + 0.5 * Fdr / N;

        //Pc
        x[2] = rho * temperature - virial / (3 * V);
        //Pm
        x[3] = 3.0 * (N - 1.0) / V * temperature * gV - virial / (3 * V) + mV * Fdr / V;
//        System.out.println(x[2] + "  " + x[3]);

        //FD
        double dbeta = 0.01;
        x[4] = 3.0 * (N - 1.0) / 2.0 / N * temperature + dbUdb(dbeta) / N;

        System.out.println(n + " " + x[0] + " " + x[1] + " " + x[4]);
//        System.out.println(n + " " + (x[1] - x[4]));

        return data;
    }


    public double dbUdb(double dbeta) {
        double beta = 1 / temperature;
        double beta_p = beta + dbeta;
        double beta_m = beta - dbeta;

        Vector dri = box.getSpace().makeVector();

        for (IAtom atom : box.getLeafList()) {
            int iAtom = atom.getLeafIndex();
            dri.Ev1Mv2(atom.getPosition(), latticePositions[iAtom]);
            box.getBoundary().nearestImage(dri);
            dri.TE(Math.sqrt(beta / beta_p));
            atom.getPosition().Ev1Pv2(latticePositions[atom.getLeafIndex()], dri);
        }
        potentialMaster.init();
        double U_p = potentialMaster.computeAll(false);

        for (IAtom atom : box.getLeafList()) {
            int iAtom = atom.getLeafIndex();
            atom.getPosition().E(rOrig[iAtom]);
        }
        potentialMaster.init();

        for (IAtom atom : box.getLeafList()) {
            int iAtom = atom.getLeafIndex();
            dri.Ev1Mv2(atom.getPosition(), latticePositions[iAtom]);
            box.getBoundary().nearestImage(dri);
            dri.TE(Math.sqrt(beta / beta_m));
            atom.getPosition().Ev1Pv2(latticePositions[atom.getLeafIndex()], dri);
        }
        potentialMaster.init();
        double U_m = potentialMaster.computeAll(false);

        for (IAtom atom : box.getLeafList()) {
            int iAtom = atom.getLeafIndex();
            atom.getPosition().E(rOrig[iAtom]);
        }
        potentialMaster.init();

        double dbUdb = (beta_p * U_p - beta_m * U_m) / 2.0 / dbeta;

        return dbUdb;
    }

    public DataTag getTag() {
        return tag;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

}