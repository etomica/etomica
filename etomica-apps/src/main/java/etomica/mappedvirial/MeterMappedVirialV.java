package etomica.mappedvirial;

 import etomica.api.IAtom;
import etomica.api.IBox;
import etomica.api.IPotentialMaster;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.iterator.IteratorDirective;
import etomica.data.DataSourceScalar;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.integrator.IntegratorVelocityVerlet.MyAgent;
import etomica.potential.PotentialCalculationForceSum;
import etomica.space.ISpace;
import etomica.units.Pressure;

public class MeterMappedVirialV extends DataSourceScalar implements  AgentSource<IntegratorVelocityVerlet.MyAgent> {

    protected final ISpace space;
    protected final IPotentialMaster potentialMaster;
    protected final PotentialCalculationForceSum pcForce;
    protected final IBox box;
    protected final IteratorDirective allAtoms;
    protected final AtomLeafAgentManager<MyAgent> forceManager;
    protected final PotentialCalculationMappedVirialV pc;
    
    public MeterMappedVirialV(ISpace space, IPotentialMaster potentialMaster, IBox box, int nbins) {
        super("pma",Pressure.DIMENSION);
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
        pc = new PotentialCalculationMappedVirialV(space, box, nbins, forceManager);
        allAtoms = new IteratorDirective();
    }
    
    public MyAgent makeAgent(IAtom a, IBox agentBox) {
        return new MyAgent(space);
    }
    
    public void releaseAgent(MyAgent agent, IAtom atom, IBox agentBox) {}

    public PotentialCalculationMappedVirialV getPotentialCalculation() {
        return pc;
    }

    public double getDataAsScalar() {
        pcForce.reset();
        potentialMaster.calculate(box, allAtoms, pcForce);
        pc.reset();
        potentialMaster.calculate(box, allAtoms, pc);
        return pc.getPressure();
    }
}
