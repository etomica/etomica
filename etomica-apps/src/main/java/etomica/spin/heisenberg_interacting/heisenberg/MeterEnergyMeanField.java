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
import etomica.potential.PotentialCalculationEnergySum;
import etomica.potential.PotentialCalculationTorqueSum;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space1d.Vector1D;
import etomica.space2d.Vector2D;
import etomica.units.dimensions.Null;
import etomica.util.numerical.BesselFunction;

public class MeterEnergyMeanField implements IDataSource, AtomLeafAgentManager.AgentSource<MeterEnergyMeanField.ForceTorque> {

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
    protected final PotentialCalculationEnergySum pce;

    public MeterEnergyMeanField(Space space, Box box, double J, PotentialMaster potentialMaster, double temperature) {
        this.potentialMaster = potentialMaster;
        this.temperature = temperature;
        pc = new PotentialCalculationMeanField(space, J, box);
        pce = new PotentialCalculationEnergySum();
        torqueSum = new PotentialCalculationTorqueSum();
        torqueAgentManager = new AtomLeafAgentManager<ForceTorque>(this, box, ForceTorque.class);
        torqueSum.setAgentManager(torqueAgentManager);
        data = new DataDoubleArray(1);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{1});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.box = box;
        allAtoms = new IteratorDirective();
    }

    public static double getpVelocity(double dtheta, double b) {
        double aC0 = FunctionCosIntegral.cosInt(dtheta, b) * BesselFunction.I(1, b) / BesselFunction.I(0, b);
        double C1 = FunctionCosIntegral.coscosInt(dtheta, b);
        return aC0 - C1;
    }


    @Override
    public IData getData() {
        pc.reset();
        torqueSum.reset();
        pce.zeroSum();

        potentialMaster.calculate(box, allAtoms, pc);
//        potentialMaster.calculate(box, allAtoms, torqueSum);
//        potentialMaster.calculate(box, allAtoms, pce);
//        double u = pce.getSum();

        AtomLeafAgentManager<Vector> agentManager = pc.getAgentManager();
        IAtomList atoms = box.getLeafList();
        data.E(0);
        double sum = 0;
        for (int i = 0; i < atoms.getAtomCount(); i++) {
            IAtomOriented a = (IAtomOriented) atoms.getAtom(i);
            Vector h = agentManager.getAgent(a);
            double hmag = Math.sqrt(h.squared());
            h.TE(1.0 / hmag);
            double I1I0 = BesselFunction.I(1, hmag / temperature) / BesselFunction.I(0, hmag / temperature);
            Vector o = a.getOrientation().getDirection();
            double sindtheta = h.getX(0) * o.getX(1) - o.getX(0) * h.getX(1);
            double cosdtheta = h.dot(o);
            double dtheta = Math.atan2(sindtheta, cosdtheta);
            double pInv = Math.exp(-hmag / temperature * cosdtheta);
            double v = getpVelocity(dtheta, hmag / temperature) * pInv;

//            sum += hmag  * (-x + cosdtheta - v * (f + hmag * sindtheta) / temperature);
//            sum += hmag  * (-x - v * (f + hmag * sindtheta) / temperature);
            sum += hmag  * (-I1I0 + v *  hmag * sindtheta / temperature);
//            usum+= hmag * cosdtheta;
        }

        data.E(sum);
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
}
