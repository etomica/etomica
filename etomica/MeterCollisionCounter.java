package etomica;
import etomica.units.*;

/**
 * This is a meter to count the number of collisions which can be used to determine efficiency vs
 * time.  This only works for a hard potential and with IntegratorHard as there is no CollisionListener
 * in the other integrators.  Average and error are useless quantities since this is merely a counter,
 * so they return NaN.
 *
 * @author Rob Riggleman
 */

public class MeterCollisionCounter extends MeterScalar implements IntegratorHardAbstract.CollisionListener {
    private int counter;
    private IntegratorHard integratorHard;
    
    public MeterCollisionCounter() {
        this(Simulation.instance);
    }
    public MeterCollisionCounter(Simulation sim) {
        super(sim);
        setLabel("Number of Collisions");
    }

    /**
     * Resets counter to zero
     */
    public void reset() {counter = 0;}
    
    public double average() {return Double.NaN;}
    
    public double error() {return Double.NaN;}
        
    public Unit defaultIOUnit() {return Count.UNIT;}
    
    public Dimension getDimension() {return Dimension.QUANTITY;}
    
    public double currentValue() {
        if(counter > 50000) {
            integratorHard.halt();
            System.out.println("Halting execution from MeterCollisionCounter");
        }
        return (double)counter;
    }
    
    
   /**
    * Implements CollisionListener.  Adds one to the counter with each collision.
    */
    public void collisionAction(IntegratorHardAbstract.Agent agent) {
        counter++;
    }
    
    /**
     * Registers itself with the integrator as CollisionListener.
     */
     
    public void setPhaseIntegrator(Integrator newIntegrator) {
        super.setPhaseIntegrator(newIntegrator);
        if(newIntegrator instanceof IntegratorHard) {
            integratorHard = (IntegratorHard)newIntegrator;
            integratorHard.addCollisionListener(this);
        }
        else throw new IllegalArgumentException("Error in integrator type in MeterCollisionCounter");
    }
}
            
    
    