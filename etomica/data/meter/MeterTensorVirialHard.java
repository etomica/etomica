package etomica.data.meter;
import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.Phase;
import etomica.Space;
import etomica.data.DataSourceCountTime;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorMD;
import etomica.space.Tensor;
import etomica.units.Dimension;

/**
 * A MeterTensor that returns the virial component of the pressure tensor for a hard potential.  
 * This is the counterpart to MeterVelocityTensor, also called by the surface tension meter.
 *
 * @author Rob Riggleman
 */

public class MeterTensorVirialHard extends MeterTensor implements IntegratorHard.CollisionListener, EtomicaElement {
    
    public MeterTensorVirialHard(Space space) {
        super(space);
        setLabel("PV/NkT");
        collisionVirial = space.makeTensor();
        virialTensor = space.makeTensor();
        virialSum = new double[space.D()][space.D()];
        timer = new DataSourceCountTime();
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Virial tensor for a hard potential");
        return info;
    }    
    
    /**
     * Dimension of the measured quantity, which is energy since virial is reported as PV/N
     */
    public Dimension getDimension() {return Dimension.NULL;}
        
    /**
     * Current value of the meter, obtained by dividing sum of collision virial contributions by time elapsed since last call.
     * If elapsed-time interval is zero, returns the value reported at the last call to the method.
     */
    //XXX phase parameter is not used appropriately here
    //TODO consider how to ensure timer is advanced before this method is invoked
    public Tensor getDataAsTensor(Phase phase) {
        double elapsedTime = timer.getData()[0];
        if(elapsedTime == 0.0) {
            virialTensor.E(Double.NaN);
            return virialTensor;
        }
        int D = phase.boundary().dimensions().D();
 
        for (int i = 0; i < collisionVirial.length(); i++) {
            for (int j = 0; j < collisionVirial.length(); j++) {
              virialTensor.setComponent(i, j, virialSum[i][j]);
              virialSum[i][j] = 0.0;
            }
        }
        virialTensor.TE(-1./(integratorHard.getTemperature()*elapsedTime*(double)(D*phase.atomCount())));
        //add 1.0 to diagonal elements
        //not included because meter returns only virial contribution to pressure
//        for (int i = 0; i < collisionVirial.length(); i++) {
//            virialTensor.PE(i,i,1.0);
//        }       
        timer.reset();
        return virialTensor;
    }
    
    /**
     * Sums contribution to virial for each collision.
     */
    public void collisionAction(IntegratorHard.Agent agent) {
       collisionVirial.E(agent.collisionPotential.lastCollisionVirialTensor());
//       collisionValue(agent);            //Assign virial
        for (int i = 0; i < collisionVirial.length(); i++) {
            for (int j = 0; j < collisionVirial.length(); j++) {
                virialSum[i][j] += collisionVirial.component(i, j);
            }
        }
    }
    
    /**
     * Contribution to the virial from the most recent collision of the given pair/potential.
     */
//    public Space.Tensor collisionValue(IntegratorHard.Agent agent) {
//        collisionVirial.E(agent.collisionPotential.lastCollisionVirialTensor());
//        collisionVirial.TE(1/(double)phase.atomCount());
//        return collisionVirial;
//    }
                
    /**
     * Informs meter of the integrator for its phase, and zeros elapsed-time counter
     */
    protected void setIntegrator(IntegratorHard newIntegrator) {
        if(newIntegrator == integratorHard) return;
        if(integratorHard != null) integratorHard.removeCollisionListener(this);
        integratorHard = newIntegrator;
        timer.setIntegrator(new IntegratorMD[] {(IntegratorMD)newIntegrator});
        if(newIntegrator != null) integratorHard.addCollisionListener(this);
    }
    
    private final DataSourceCountTime timer;
    /**
     * Sums contributions to the virial from collisions, between calls to currentValue
     */
    private final double[][] virialSum;
    /**
     * Tensor used to return current value of the meter.
     */
    private final Tensor virialTensor;
    /**
     * Tensor used to return meter's collision value (value contributed from last collision).
     */
    private Tensor collisionVirial;
    /**
     * Convenience handle to integrator to avoid multiple casting to IntegratorHard
     */
    private IntegratorHard integratorHard;

}//end of MeterTensorVirialHard