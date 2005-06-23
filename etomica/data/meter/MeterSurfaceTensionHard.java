package etomica.data.meter;
import etomica.DataInfo;
import etomica.EtomicaInfo;
import etomica.Meter;
import etomica.Phase;
import etomica.Space;
import etomica.data.DataSourceScalar;
import etomica.space.Tensor;
import etomica.units.Dimension;

/**
 * This is a meter to measure the surface tension for a hard potential.  
 * It uses the meter that measures the pressure tensor to do this
 * and returns the result as PV/N, which has dimensions of energy.
 * Assumes a slab geometry for the liquid, so surface tension is divided by 2 to account for both interfaces.
 *
 * @author Rob Riggleman
 */

public class MeterSurfaceTensionHard extends DataSourceScalar implements Meter {//, EtomicaElement {
    
    public MeterSurfaceTensionHard(Space space) {
        super(new DataInfo("Surface Tension",Dimension.ENERGY));
        velocityTensor = new MeterTensorVelocity(space);
        virialTensor = new MeterTensorVirialHard(space);
        pressureTensor = space.makeTensor();
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
    
    public double getDataAsScalar() {
        if (phase == null) throw new IllegalStateException("must call setPhase before using meter");
        pressureTensor.E(velocityTensor.getDataAsTensor(phase));
        pressureTensor.PE(virialTensor.getDataAsTensor(phase));
        switch (phase.space().D) {
            case 1:
                surfaceTension = pressureTensor.component(0, 0);
                break;
            case 2:
                surfaceTension = 0.5*(pressureTensor.component(0, 0) - pressureTensor.component(1, 1));
                break;
            case 3:
                surfaceTension = 0.5*(pressureTensor.component(0, 0) - 0.5*(pressureTensor.component(1, 1) + pressureTensor.component(2, 2)));
                break;
        }
        return surfaceTension;
    }
    
    /**
     * @return Returns the phase.
     */
    public Phase getPhase() {
        return phase;
    }
    /**
     * @param phase The phase to set.
     */
    public void setPhase(Phase phase) {
        if (phase.space().D < 1 || phase.space().D > 3) throw new IllegalArgumentException("phase's space must be 1, 2 or 3 dimensional");
        this.phase = phase;
    }

    private Phase phase;
    private final Tensor pressureTensor;
    private final MeterTensorVelocity velocityTensor;
    private final MeterTensorVirialHard virialTensor;
    private double surfaceTension;
}
