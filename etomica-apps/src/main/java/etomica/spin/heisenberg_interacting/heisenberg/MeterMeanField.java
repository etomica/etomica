package etomica.spin.heisenberg_interacting.heisenberg;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
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

    public MeterMeanField(Space space, Box box, double J, PotentialMaster potentialMaster, double temperature) {
        this.potentialMaster = potentialMaster;
        this.temperature = temperature;
        pc = new PotentialCalculationMeanField(space, J, box);
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

    @Override
    public IData getData() {
        pc.reset();
        torqueSum.reset();

        potentialMaster.calculate(box, allAtoms, pc);
        potentialMaster.calculate(box, allAtoms, torqueSum);

        AtomLeafAgentManager<Vector> agentManager = pc.getAgentManager();
        IAtomList atoms = box.getLeafList();
        for (int i = 0; i < atoms.getAtomCount(); i++) {
            IAtom a = atoms.getAtom(i);
            Vector h = agentManager.getAgent(a);
            double hmag = Math.sqrt(h.squared());
            double x = BesselFunction.doCalc(true, 1, hmag / temperature) / BesselFunction.doCalc(true, 0, hmag / temperature);
            h.TE(x / hmag);
            double f = torqueAgentManager.getAgent(a).torque().getX(0);

        }
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
