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
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Null;
import etomica.util.numerical.BesselFunction;

public class MeterMeanField implements IDataSource {

    protected final DataDoubleArray data;
    protected final DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected final PotentialCalculationMeanField pc;
    protected final PotentialMaster potentialMaster;
    protected final Box box;
    protected final IteratorDirective allAtoms;
    protected double temperature;

    public MeterMeanField(Space space, Box box, double J, PotentialMaster potentialMaster, double temperature) {
        this.potentialMaster = potentialMaster;
        this.temperature = temperature;
        pc = new PotentialCalculationMeanField(space, J, box);
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
        potentialMaster.calculate(box, allAtoms, pc);
        AtomLeafAgentManager<Vector> agentManager = pc.getAgentManager();
        IAtomList atoms = box.getLeafList();
        for (int i = 0; i < atoms.getAtomCount(); i++) {
            IAtom a = atoms.getAtom(i);
            Vector h = agentManager.getAgent(a);
            double hmag = Math.sqrt(h.squared());
            double x = BesselFunction.doCalc(true, 1, hmag / temperature) / BesselFunction.doCalc(true, 0, hmag / temperature);
            h.TE(x / hmag);


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
}
