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
import etomica.space.Vector;
import etomica.space1d.Vector1D;
import etomica.space2d.Vector2D;
import etomica.units.dimensions.Null;

public class MeterDipoleMoment implements IDataSource, AtomLeafAgentManager.AgentSource<MeterDipoleMoment.ForceTorque> {

    protected final DataDoubleArray data;
    protected final DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected final Box box;

    public MeterDipoleMoment(Box box) {
        data = new DataDoubleArray(2);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{2});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.box = box;
    }

    @Override
    public IData getData() {
        IAtomList atoms = box.getLeafList();
        data.E(0);
        double[] d = data.getData();
        for (int i = 0; i < atoms.size(); i++) {
            IAtomOriented a = (IAtomOriented) atoms.get(i);
            Vector o = a.getOrientation().getDirection();

            d[0] += o.getX(0);
            d[1] += o.getX(1);
        }
        data.TE(1.0 / atoms.size());
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
