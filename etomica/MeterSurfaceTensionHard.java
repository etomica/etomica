package simulate;
import simulate.units.*;

/**
 * This is a meter to measure the surface tension for a hard potential.  
 * It uses the meter that measures the pressure tensor to do this
 * and returns the result as PV/N, which has dimensions of energy.
 * Assumes a slab geometry for the liquid, so surface tension is divided by 2 to account for both interfaces.
 *
 * @author Rob Riggleman
 */

public class MeterSurfaceTensionHard extends Meter implements Meter.Atomic, Meter.Collisional {
    private final Space.Tensor pressureTensor;
    private final MeterTensorVelocity velocityTensor = new MeterTensorVelocity();
    private final MeterTensorVirialHard virialTensor = new MeterTensorVirialHard();
    private double surfaceTension, collisionValue, velocityValue;
    private int D;
    
    public MeterSurfaceTensionHard() {
        this(Simulation.instance);
    }
    public MeterSurfaceTensionHard(Simulation sim) {
        super(sim);
        pressureTensor = sim.space().makeTensor();
        D = sim.space().D();
        setLabel("Surface Tension");
    }
    
    /**
     * Returns dimensions of this meters measured quanitity, which in this case is energy
     */
    public final Dimension getDimension() {return Dimension.ENERGY;}
    
    /**
     * Indicates that this meter does not reference the phase boundary.
     * @return false
     */
    public final boolean usesPhaseBoundary() {return false;}
    
    /**
     * Indicates that this meter does not use any iterators.
     * @return false
     */
    public final boolean usesPhaseIteratorFactory() {return false;}
    
    /**
     * Gives current value of the surface tension, obtained by summing velocity and virial contributions.
     * Virial contribution includes sum over all collisions since last call to this method, while
     * velocity contribution is based on atom velocities in current configuration.
     * Surface tension is given by difference between normal and tangential stress components.
     * Written for 2-, or 3-dimensional systems; assumes that normal to interface is along x-axis.
     */
    public final double currentValue() {
        pressureTensor.E(velocityTensor.currentValue());
        pressureTensor.PE(virialTensor.currentValue());
        if (D == 1) {surfaceTension = pressureTensor.component(0, 0);}
        else if (D == 2) {surfaceTension = 0.5*(pressureTensor.component(0, 0) - pressureTensor.component(1, 1));}
        else {surfaceTension = 0.5*(pressureTensor.component(0, 0) - 0.5*(pressureTensor.component(1, 1) + pressureTensor.component(2, 2)));}
        return surfaceTension;
    }
    
    /**
     * Calls collisionAction method of virial contribution
     */
    public void collisionAction(AtomPair pair, Potential.Hard p) {
        virialTensor.collisionAction(pair, p);
    }
    
    /**
     * Contribution to the surface tension from the recent collision of the given pair
     */
    public double collisionValue(AtomPair pair, Potential.Hard p) {
        Space.Tensor vTensor = virialTensor.collisionValue(pair, p);
        if (D == 1) {collisionValue = vTensor.component(0, 0);}
        else if (D == 2) {collisionValue = vTensor.component(0, 0) - vTensor.component(1, 1);}
        else {collisionValue = vTensor.component(0, 0) - 0.5*(vTensor.component(1, 1) + vTensor.component(2, 2));}
        return 0.5*collisionValue;
    }
    
    /**
     * Calls current value method of the velocity tensor
     */
    public final double currentValue(Atom a) {
        Space.Tensor vTensor = velocityTensor.currentValue(a);
        if (D == 1) {velocityValue = vTensor.component(0, 0);}
        else if (D == 2) {velocityValue = vTensor.component(0, 0) - vTensor.component(1, 1);}
        else {velocityValue = vTensor.component(0, 0) - 0.5*(vTensor.component(1, 1) + vTensor.component(2, 2));}
        return 0.5*velocityValue;
    }
}