package etomica;
import etomica.units.*;

public class MeterPressureHardTensor extends MeterTensor implements IntegratorHardAbstract.CollisionListener {
    
    private double[][] virialSum;
    private double t0, t, velocity2;
    private final AtomIteratorList ai1 = new AtomIteratorList();
    private Space.Tensor velocityTensor, v, pressureTensor;
    private IntegratorHard integratorHard;
    private int D;
    
    public MeterPressureHardTensor() {
        this(Simulation.instance);
    }
    public MeterPressureHardTensor(Simulation sim) {
        super(sim);
        D = sim.space().D();
        virialSum = new double[D][D];
        pressureTensor = sim.space().makeTensor();
        velocityTensor = sim.space().makeTensor();
        v = sim.space().makeTensor();
        setLabel("PV/Nk");
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Pressure tensor measured via components of virial averaged over hard collisions");
        return info;
    }
    
    public void setPhase(Phase p) {
        super.setPhase(p);
        ai1.setBasis(p.speciesMaster.atomList);
    }

    public final boolean usesPhaseBoundary() {return false;}
     
    public Dimension getDimension() {return Dimension.TEMPERATURE;}
    
    public Space.Tensor currentValue() {
        double t = integratorHard.elapsedTime();
        if (t > t0) {
            velocityTensor.E(0.);
            ai1.reset();
            while (ai1.hasNext()) {
                Atom a = ai1.next();
                v.E(a.coord.momentum(), a.coord.momentum());
                v.TE((a.coord.rm()));
                velocityTensor.PE(v);
            }
            
            velocityTensor.TE(1.0/phase.atomCount());
                    
            for (int h=0; h<D; h++) {
                for (int j=0; j<D; j++) {
                    pressureTensor.setComponent(h, j, virialSum[h][j]);
                    virialSum[h][j] = 0.0;
                }
            }
        
            pressureTensor.TE(-1/((t-t0)*(double)(D*phase.atomCount())));
            pressureTensor.PE(velocityTensor);
        
            t0 = t;
        }
        return pressureTensor;
    }
    
    public double currentValue(int i, int j) {
        return currentValue().component(i, j);
    }
    
    public void collisionAction(IntegratorHardAbstract.Agent agent) {
        Space.Tensor lcvt = agent.collisionPotential.lastCollisionVirialTensor();
        for (int i=0; i<lcvt.length(); i++) {
            for (int j=0; j<lcvt.length(); j++) {
                virialSum[i][j] += lcvt.component(i, j);
            }
        }
    }
    
    protected void setPhaseIntegrator(Integrator newIntegrator) {
        super.setPhaseIntegrator(newIntegrator);
        if(newIntegrator instanceof IntegratorHard) {
            integratorHard = (IntegratorHard)newIntegrator;
            integratorHard.addCollisionListener(this);
            t0 = integratorHard.elapsedTime();
        }
        else {
            System.out.println("Error in integrator type in MeterPressureHardTensor");
            System.exit(1);
        }
    }
}