package etomica.spin.heisenberg_interacting.heisenberg;

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
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationTorqueSum;
import etomica.potential.PotentialMaster;
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
    protected final PotentialCalculationTorqueSum torqueSum;
    protected final PotentialCalculationMeanField pc;
    protected final PotentialMaster potentialMaster;
    protected final Box box;
    protected final IteratorDirective allAtoms;
    protected double temperature;
    protected final AtomLeafAgentManager<ForceTorque> torqueAgentManager;
    protected final List<Vector> spins;

    public MeterMeanField(Space space, Box box, double J, PotentialMaster potentialMaster, double temperature) {
        this.potentialMaster = potentialMaster;
        this.temperature = temperature;
        pc = new PotentialCalculationMeanField(space, J, box);
        torqueSum = new PotentialCalculationTorqueSum();
        torqueAgentManager = new AtomLeafAgentManager<ForceTorque>(this, box, ForceTorque.class);
        torqueSum.setAgentManager(torqueAgentManager);
        data = new DataDoubleArray(2);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{2});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.box = box;
        allAtoms = new IteratorDirective();
        spins = new ArrayList<>();
    }

    public static double[] getVelocity(double dtheta, double b, double cosTheta0, double sinTheta0) {
        double aC0 = FunctionCosIntegral.cosInt(dtheta, b) * BesselFunction.I(1, b) / BesselFunction.I(0, b);
        double C1 = FunctionCosIntegral.coscosInt(dtheta, b);
        double S1 = FunctionCosIntegral.sincosInt(dtheta, b);
        double A = Math.exp(-b * Math.cos(dtheta));
        return new double[]{A * (aC0 * cosTheta0 - C1 * cosTheta0 + S1 * sinTheta0),
                A * (aC0 * sinTheta0 - C1 * sinTheta0 + S1 * cosTheta0)};
    }


    @Override
    public IData getData() {
        pc.reset();
        torqueSum.reset();

        potentialMaster.calculate(box, allAtoms, pc);
        potentialMaster.calculate(box, allAtoms, torqueSum);

        AtomLeafAgentManager<Vector> agentManager = pc.getAgentManager();
        IAtomList atoms = box.getLeafList();
        if (spins.size() < atoms.getAtomCount()) {
            for (int i = spins.size(); i < atoms.getAtomCount(); i++) {
                spins.add(new Vector2D());
            }
        }
        data.E(0);
        double[] d = data.getData();
        for (int i = 0; i < atoms.getAtomCount(); i++) {
            IAtomOriented a = (IAtomOriented) atoms.getAtom(i);
            Vector h = agentManager.getAgent(a);
            double hmag = Math.sqrt(h.squared());
            h.TE(1.0 / hmag);
            double x = BesselFunction.doCalc(true, 1, hmag / temperature) / BesselFunction.doCalc(true, 0, hmag / temperature);
            double f = torqueAgentManager.getAgent(a).torque().getX(0);
            Vector o = a.getOrientation().getDirection();
            double sindtheta = h.getX(0) * o.getX(1) - o.getX(0) * h.getX(1);
            double cosdtheta = h.dot(o);
            double y = (f + hmag * sindtheta) / temperature;
            double dtheta = Math.atan2(sindtheta, cosdtheta);
            double theta0 = Math.atan2(h.getX(1), h.getX(0));
            double[] v = getVelocity(dtheta, hmag / temperature, h.getX(0), h.getX(1));
            double[] vc = getVelocity(0 - theta0, hmag / temperature, h.getX(0), h.getX(1));
            double[] vs = getVelocity(Math.PI / 2 - theta0, hmag / temperature, h.getX(0), h.getX(1));
            Vector s = spins.get(i);
            s.setX(0, h.getX(0) * x + y * (v[0] - vc[0]));
            s.setX(1, h.getX(1) * x + y * (v[1] - vs[1]));

            d[0] += s.getX(0);
            d[1] += s.getX(1);
        }
        data.TE(1.0 / atoms.getAtomCount());
        return data;
    }

    public Vector getSpin(IAtom a) {
        return spins.get(a.getLeafIndex());
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
