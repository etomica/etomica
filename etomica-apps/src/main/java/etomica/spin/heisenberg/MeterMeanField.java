package etomica.spin.heisenberg;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.integrator.Integrator;
import etomica.potential.compute.NeighborIterator;
import etomica.potential.compute.NeighborManager;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space1d.Vector1D;
import etomica.space2d.Vector2D;
import etomica.units.dimensions.Null;
import etomica.util.numerical.BesselFunction;

import java.util.ArrayList;
import java.util.List;

public class MeterMeanField implements IDataSource, AtomLeafAgentManager.AgentSource<MeterMeanField.ForceTorque> {

    protected final DataDoubleArray data;
    protected final DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected final NeighborIterator nbrIterator;
    protected final Box box;
    protected double temperature;
    protected double[] I1I0, I2etc, eta, cost0, sint0, cosdtheta, sindtheta, dtheta, theta0;
    protected final List<Vector> spins, thetaDot;
    protected final double J2;
    protected final Vector[] netOrientation;

    public MeterMeanField(Space space, Box box, double J, NeighborManager nbrManager, double temperature) {
        this.nbrIterator = nbrManager.makeNeighborIterator();
        this.temperature = temperature;
        this.J2 = 0.5*J;
        data = new DataDoubleArray(2);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{2});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.box = box;
        spins = new ArrayList<>();
        thetaDot = new ArrayList<>();
        netOrientation = new Vector[box.getLeafList().size()];
        for (int i=0; i<netOrientation.length; i++) {
            netOrientation[i] = box.getSpace().makeVector();
        }
    }

    public static double[] getpVelocity(double dtheta, double b, double cosTheta0, double sinTheta0) {
        double aC0 = FunctionCosIntegral.cosInt(dtheta, b) * BesselFunction.I(1, b) / BesselFunction.I(0, b);
        double C1 = FunctionCosIntegral.coscosInt(dtheta, b);
        double S1 = FunctionCosIntegral.sincosInt(dtheta, b);
        return new double[]{(aC0 * cosTheta0 - C1 * cosTheta0 + S1 * sinTheta0),
                (aC0 * sinTheta0 - C1 * sinTheta0 - S1 * cosTheta0)};
    }


    @Override
    public IData getData() {
        for (int i=0; i<netOrientation.length; i++) {
            netOrientation[i].E(0);
        }
        for (IAtom a1 : box.getLeafList()) {
            IAtomOriented atom1 = (IAtomOriented) a1;
            nbrIterator.iterUpNeighbors(a1.getLeafIndex(), new NeighborIterator.NeighborConsumer() {
                @Override
                public void accept(IAtom jAtom, Vector rij, int n) {
                    IAtomOriented atom2 = (IAtomOriented) jAtom;
                    Vector ei = atom1.getOrientation().getDirection();
                    Vector ej = atom2.getOrientation().getDirection();

                    netOrientation[atom1.getLeafIndex()].PEa1Tv1(J2, ej);
                    netOrientation[atom2.getLeafIndex()].PEa1Tv1(J2, ei);
                }
            });
        }

        IAtomList atoms = box.getLeafList();
        int atomCount = atoms.size();
        if (spins.size() < atomCount) {
            for (int i = spins.size(); i < atomCount; i++) {
                spins.add(new Vector2D());
                thetaDot.add(new Vector2D());
            }
            I1I0 = new double[atomCount];
            I2etc = new double[atomCount];
            eta = new double[atomCount];
            cost0 = new double[atomCount];
            sint0 = new double[atomCount];
            cosdtheta = new double[atomCount];
            sindtheta = new double[atomCount];
            dtheta = new double[atomCount];
            theta0 = new double[atomCount];
        }
        data.E(0);
        double[] d = data.getData();
        for (int i = 0; i < atoms.size(); i++) {
            IAtomOriented a = (IAtomOriented) atoms.get(i);
            Vector h = netOrientation[i];
            eta[i] = Math.sqrt(h.squared());
            cost0[i] = h.getX(0) / eta[i];
            sint0[i] = h.getX(1) / eta[i];
            double bh = eta[i] / temperature;

            double I0 = BesselFunction.I(0, bh);
            double I1 = BesselFunction.I(1, bh);
            double I2 = BesselFunction.I(2, bh);
            I1I0[i] = I1 / I0;
            I2etc[i] = 0.5 * (I2 / I0 - 2 * I1I0[i] * I1I0[i] + 1);

            Vector o = a.getOrientation().getDirection();
            sindtheta[i] = cost0[i] * o.getX(1) - sint0[i] * o.getX(0);
            cosdtheta[i] = cost0[i] * o.getX(0) + sint0[i] * o.getX(1);
            dtheta[i] = Math.atan2(sindtheta[i], cosdtheta[i]);
            theta0[i] = Math.atan2(sint0[i], cost0[i]);
            double[] v = getpVelocity(dtheta[i], bh, cost0[i], sint0[i]);
            double[] vc = getpVelocity(0 - theta0[i], bh, cost0[i], sint0[i]);
            double[] vs = getpVelocity(Math.PI / 2 - theta0[i], bh, cost0[i], sint0[i]);
            double pInv = Math.exp(-bh * cosdtheta[i]);
            Vector s = spins.get(i);
            Vector tDot = thetaDot.get(i);
            tDot.setX(0,pInv * (v[0] - vc[0]));
            tDot.setX(1,pInv * (v[1] - vs[1]));
            s.setX(0, cost0[i] * I1I0[i] - bh * sindtheta[i] * tDot.getX(0));//mapped-average estimate of spin orientation
            s.setX(1, sint0[i] * I1I0[i] - bh * sindtheta[i] * tDot.getX(1));

            d[0] += s.getX(0);
            d[1] += s.getX(1);
        }
        data.TE(1.0 / atoms.size());
        return data;
    }

    public Vector getSpin(IAtom a) {
        return spins.get(a.getLeafIndex());
    }

    public Vector getThetaDot(IAtom a) {return thetaDot.get(a.getLeafIndex()); }

    public double[] getI1I0() {return I1I0;}

    public double[] getI2etc() {return I2etc;}

    public double[] getEta() {return eta;}

    public double[] getCost0() {return cost0;}

    public double[] getSint0() {return sint0;}

    public double[] getSindtheta() {return sindtheta;}

    public double[] getCosdtheta() {return cosdtheta;}

    public double[] getDtheta() {return dtheta;}

    public double[] getTheta0() {return theta0;}


    @Override
    public DataTag getTag() {
        return tag;
    }

    @Override
    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    @Override
    public ForceTorque makeAgent(IAtom a, Box agentBox) {
        return new ForceTorque();
    }

    @Override
    public void releaseAgent(ForceTorque agent, IAtom atom, Box agentBox) {

    }

    public static class ForceTorque implements Integrator.Forcible, Integrator.Torquable {
        protected final Vector f, t;

        public ForceTorque() {
            f = new Vector2D();
            t = new Vector1D();
        }

        public Vector force() {
            return f;
        }

        public Vector torque() {
            return t;
        }
    }
}
