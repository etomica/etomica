package etomica;
import etomica.units.Dimension;

/**
 * This is a meter to measure the surface tension for a hard potential.  
 * It uses the meter that measures the pressure tensor to do this
 * and returns the result as PV/N, which has dimensions of energy.
 * Assumes a slab geometry for the liquid, so surface tension is divided by 2 to account for both interfaces.
 *
 * @author Rob Riggleman
 */

public class MeterSurfaceTensionHard extends MeterScalar {//, EtomicaElement {
    private final Space.Tensor pressureTensor;
    private final MeterTensorVelocity velocityTensor;
    private final MeterTensorVirialHard virialTensor;
    private double surfaceTension, collisionValue, velocityValue;
    private int D;
    
    public MeterSurfaceTensionHard(Space space) {
        super();
        velocityTensor = new MeterTensorVelocity(space);
        virialTensor = new MeterTensorVirialHard(space);
        pressureTensor = space.makeTensor();
        D = space.D();
        setLabel("Surface Tension");
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Measures surface tension for a collision-based potential");
        return info;
    }

    /**
     * Returns dimensions of this meters measured quanitity, which in this case is energy
     */
    public final Dimension getDimension() {return Dimension.ENERGY;}
    
    /**
     * Gives current value of the surface tension, obtained by summing velocity and virial contributions.
     * Virial contribution includes sum over all collisions since last call to this method, while
     * velocity contribution is based on atom velocities in current configuration.
     * Surface tension is given by difference between normal and tangential stress components.
     * Written for 2-, or 3-dimensional systems; assumes that normal to interface is along x-axis.
     */
    
    public double getDataAsScalar(Phase p) {
        pressureTensor.E(velocityTensor.getDataAsTensor(p));
        pressureTensor.PE(virialTensor.getDataAsTensor(p));
        if (D == 1) {surfaceTension = pressureTensor.component(0, 0);}
        else if (D == 2) {surfaceTension = 0.5*(pressureTensor.component(0, 0) - pressureTensor.component(1, 1));}
        else {surfaceTension = 0.5*(pressureTensor.component(0, 0) - 0.5*(pressureTensor.component(1, 1) + pressureTensor.component(2, 2)));}
        return surfaceTension;
    }
    
}
