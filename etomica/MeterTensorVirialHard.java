package etomica;
import etomica.units.*;

/**
 * A MeterTensor that returns the virial component of the pressure tensor for a hard potential.  
 * This is the counterpart to MeterVelocityTensor, also called by the surface tension meter.
 *
 * @author Rob Riggleman
 */

/* History
 * 03/12/04 (DAK) modified collisionValue method to not check explicitly
 * for instance of PotentialNull
 */

public class MeterTensorVirialHard extends MeterTensor implements IntegratorHard.CollisionListener, EtomicaElement {
    
    /**
     * Sums contributions to the virial from collisions, between calls to currentValue
     */
    private final double[][] virialSum;
    /**
     * Tensor used to return current value of the meter.
     */
    private final Space.Tensor virialTensor;
    /**
     * Tensor used to return meter's collision value (value contributed from last collision).
     */
    private Space.Tensor collisionVirial;
    /**
     * Convenience handle to integrator to avoid multiple casting to IntegratorHard
     */
    private IntegratorHard integratorHard;
    /**
     * Simulation time of last call to currentValue.  Used to determine elapsed time for summing of collision virials.
     */
    private double t0;
    
    public MeterTensorVirialHard() {
        this(Simulation.instance);
    }
    public MeterTensorVirialHard(Simulation sim) {
        super(sim);
        collisionVirial = sim.space().makeTensor();
        virialTensor = sim.space().makeTensor();
        virialSum = new double[sim.space().D()][sim.space().D()];
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Virial tensor for a hard potential");
        return info;
    }    
    
    /**
     * Dimension of the measured quantity, which is energy since virial is reported as PV/N
     */
    public Dimension getDimension() {return Dimension.ENERGY;}
    
    /**
     * Descriptive label
     */
    public String getLabel() {return "PV/N";}
    
    /**
     * Current value of the meter, obtained by dividing sum of collision virial contributions by time elapsed since last call.
     * If elapsed-time interval is zero, returns the value reported at the last call to the method.
     */
    public Space.Tensor getData() {
        double t = integratorHard.elapsedTime();
        if (t > t0) {
            for (int i = 0; i < collisionVirial.length(); i++) {
                for (int j = 0; j < collisionVirial.length(); j++) {
                  virialTensor.setComponent(i, j, virialSum[i][j]);
                  virialSum[i][j] = 0.0;
                }
            }
            virialTensor.TE(-1./((t-t0)));
        }
        
        t0 = t;
        return virialTensor;
    }
    
    /**
     * Sums contribution to virial for each collision.
     */
    public void collisionAction(IntegratorHard.Agent agent) {
        collisionValue(agent);            //Assign virial
        for (int i = 0; i < collisionVirial.length(); i++) {
            for (int j = 0; j < collisionVirial.length(); j++) {
                virialSum[i][j] += collisionVirial.component(i, j);
            }
        }
    }
    
    /**
     * Contribution to the virial from the most recent collision of the given pair/potential.
     */
    public Space.Tensor collisionValue(IntegratorHard.Agent agent) {
        collisionVirial.E(agent.collisionPotential.lastCollisionVirialTensor());
        collisionVirial.TE(1/(double)phase.atomCount());
       return collisionVirial;
    }
                
    /**
     * Informs meter of the integrator for its phase, and zeros elapsed-time counter
     */
    protected void setPhaseIntegrator(Integrator newIntegrator) {
        super.setPhaseIntegrator(newIntegrator);
        if(newIntegrator instanceof IntegratorHard) {
            integratorHard = (IntegratorHard)newIntegrator;
            integratorHard.addCollisionListener(this);
            t0 = integratorHard.elapsedTime();
        }
        else throw new IllegalArgumentException("Error in integrator type in MeterPressureHardTensor");
    }
}//end of MeterTensorVirialHard