package etomica.mappedvirial;

 import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.potential.PotentialMaster;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.iterator.IteratorDirective;
 import etomica.data.DataSourceScalar;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.integrator.IntegratorVelocityVerlet.MyAgent;
import etomica.potential.PotentialCalculationForceSum;
 import etomica.space.Space;
 import etomica.units.Pressure;

public class MeterMappedVirialV extends DataSourceScalar implements  AgentSource<IntegratorVelocityVerlet.MyAgent> {

    protected final Space space;
    protected final PotentialMaster potentialMaster;
    protected final PotentialCalculationForceSum pcForce;
    protected final Box box;
    protected final IteratorDirective allAtoms;
    protected final AtomLeafAgentManager<MyAgent> forceManager;
    protected final PotentialCalculationMappedVirialV pc;
    
    public MeterMappedVirialV(Space space, PotentialMaster potentialMaster, Box box, int nbins) {
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
    
    public MyAgent makeAgent(IAtom a, Box agentBox) {
        return new MyAgent(space);
    }
    
    public void releaseAgent(MyAgent agent, IAtom atom, Box agentBox) {}

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
