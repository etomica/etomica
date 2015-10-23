package etomica.mappedvirial;

 import etomica.api.IAtom;
import etomica.api.IBox;
import etomica.api.IPotentialMaster;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.iterator.IteratorDirective;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.integrator.IntegratorVelocityVerlet.MyAgent;
import etomica.potential.PotentialCalculationForceSum;
import etomica.space.ISpace;
import etomica.units.Pressure;

public class MeterMappedVirial implements IEtomicaDataSource, AgentSource<IntegratorVelocityVerlet.MyAgent> {

    protected final ISpace space;
    protected final IPotentialMaster potentialMaster;
    protected final PotentialCalculationForceSum pcForce;
    protected final IBox box;
    protected final IteratorDirective allAtoms;
    protected final AtomLeafAgentManager<MyAgent> forceManager;
    protected final DataDoubleArray data;
    protected final DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected final PotentialCalculationMappedVirial pc;
    
    public MeterMappedVirial(ISpace space, IPotentialMaster potentialMaster, double pCut, IBox box, int nbins, double[] rCutoff) {
        this.space = space;
        this.box = box;
        this.potentialMaster = potentialMaster;
        pcForce = new PotentialCalculationForceSum();
        if (box != null) {
            forceManager = new AtomLeafAgentManager<MyAgent>(this, box, MyAgent.class);
            pcForce.setAgentManager(forceManager);
        }
        else {
            forceManager = null;
        }
        pc = new PotentialCalculationMappedVirial(space, box, nbins, rCutoff, forceManager, pCut);
        allAtoms = new IteratorDirective();
        data = new DataDoubleArray(rCutoff.length);
        dataInfo = new DataInfoDoubleArray("mapped virial", Pressure.DIMENSION, new int[]{rCutoff.length});
        tag = new DataTag();
        dataInfo.addTag(tag);
    }
    
    public MyAgent makeAgent(IAtom a) {
        return new MyAgent(space);
    }
    
    public void releaseAgent(MyAgent agent, IAtom atom) {}

    public PotentialCalculationMappedVirial getPotentialCalculation() {
        return pc;
    }

    public IData getData() {
        pcForce.reset();
        potentialMaster.calculate(box, allAtoms, pcForce);
        pc.reset();
        potentialMaster.calculate(box, allAtoms, pc);
        double[] pressure = pc.getPressure();
        
        double[] x = data.getData();
        for (int j=0; j<x.length; j++) {
            x[j] = pressure[j];
        }
        return data;
    }

    public DataTag getTag() {
        return tag;
    }

    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }
}
