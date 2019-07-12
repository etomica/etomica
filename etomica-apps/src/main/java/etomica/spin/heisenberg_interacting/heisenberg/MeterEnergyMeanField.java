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
import etomica.potential.*;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space1d.Vector1D;
import etomica.space2d.Vector2D;
import etomica.units.dimensions.Null;
import etomica.util.numerical.BesselFunction;

public class MeterEnergyMeanField implements IDataSource, AtomLeafAgentManager.AgentSource<MeterEnergyMeanField.ForceTorque> {

    protected final double J;
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
    protected final PotentialCalculationStuff pcExtra;
    protected double[] thetadot = new double[0];
    protected double cvSumExtra;

    public MeterEnergyMeanField(Space space, Box box, double J, PotentialMaster potentialMaster, double temperature) {
        this.potentialMaster = potentialMaster;
        this.temperature = temperature;
        this.J = J;
        pc = new PotentialCalculationMeanField(space, J, box);
        pcExtra = new PotentialCalculationStuff();
        torqueSum = new PotentialCalculationTorqueSum();
        torqueAgentManager = new AtomLeafAgentManager<ForceTorque>(this, box, ForceTorque.class);
        torqueSum.setAgentManager(torqueAgentManager);
        data = new DataDoubleArray(2);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{2});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.box = box;
        allAtoms = new IteratorDirective();
    }

    @Override
    public IData getData() {
        if (thetadot.length != box.getLeafList().getAtomCount()) {
            thetadot = new double[box.getLeafList().getAtomCount()];
        }
        pc.reset();

        potentialMaster.calculate(box, allAtoms, pc);

        AtomLeafAgentManager<Vector> agentManager = pc.getAgentManager();
        IAtomList atoms = box.getLeafList();
        data.E(0);
        double sum = 0;
        double cvSum = 0;
        for (int i = 0; i < atoms.getAtomCount(); i++) {
            IAtomOriented a = (IAtomOriented) atoms.getAtom(i);
            Vector h = agentManager.getAgent(a);
            double hmag = Math.sqrt(h.squared());
            h.TE(1.0 / hmag);
            double I1 = BesselFunction.I(1, hmag / temperature);
            double I0 = BesselFunction.I(0, hmag / temperature);
            double I2 = BesselFunction.I(2, hmag / temperature);
            Vector o = a.getOrientation().getDirection();
            double sindtheta = h.getX(0) * o.getX(1) - o.getX(0) * h.getX(1);
            double cosdtheta = h.dot(o);
            double dtheta = Math.atan2(sindtheta, cosdtheta);
            double pInv = Math.exp(-hmag / temperature * cosdtheta);

            double I1I0 = I1 / I0;
            double C0 = FunctionCosIntegral.cosInt(dtheta, hmag / temperature);
            double C1 = FunctionCosIntegral.coscosInt(dtheta, hmag / temperature);

            double v = (C0 * I1I0 - C1) * pInv;
            thetadot[i] = v;

//            sum += hmag  * (-x + cosdtheta - v * (f + hmag * sindtheta) / temperature);
//            sum += hmag  * (-x - v * (f + hmag * sindtheta) / temperature);
            double x = hmag * (-I1I0 + v * hmag * sindtheta / temperature);
            sum += x;
//            usum+= hmag * cosdtheta;

            double vbb = hmag * hmag * ((C1 * (I1I0 + cosdtheta) - C0 * I1I0 - FunctionCosIntegral.cos2cosInt(dtheta, hmag / temperature))
                    + (I2 / I0 - 2 * I1I0 * I1I0 + 1) * C0) * pInv;

            cvSum += hmag * hmag * 0.5 * (I2 / I0 - 2 * I1I0 * I1I0 + 1);
            cvSum += hmag / temperature * (-vbb + x - hmag * cosdtheta + 5 * temperature) * sindtheta
                    + v * cosdtheta - v * v * 2 * hmag * cosdtheta;
        }

        cvSumExtra = 0;
        potentialMaster.calculate(box, allAtoms, pc);
        cvSum += J * cvSumExtra;
        double[] y = data.getData();

        y[0] = sum;
        y[1] = sum * sum + cvSum;
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

    public class PotentialCalculationStuff implements PotentialCalculation {
        @Override
        public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
            IAtomOriented iatom = (IAtomOriented) atoms.getAtom(0);
            IAtomOriented jatom = (IAtomOriented) atoms.getAtom(1);
            Vector io = iatom.getOrientation().getDirection();
            Vector jo = jatom.getOrientation().getDirection();
            cvSumExtra += io.dot(jo) * thetadot[iatom.getLeafIndex()] * thetadot[jatom.getLeafIndex()];
        }
    }
}
