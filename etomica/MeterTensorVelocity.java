package etomica;
import etomica.units.*;

/**
 * A meter to compute the velocity component of the pressure tensor. 
 * Averages a tensor quantity formed from a dyad of the velocity of each atom. 
 * Specifically, the quantity averaged is 1/N * sum(pp/m), where p is the momentum,
 * m is the mass, and the sum is over all N atoms.
 * 
 * @author Rob Riggleman
 */
public class MeterTensorVelocity extends MeterTensor implements MeterTensor.Atomic {
    /**
     * Iterator of atoms.
     */
    private Atom.Iterator ai1;
    /**
     * Tensor used to form velocity dyad for each atom, and returned by currentValue(atom) method.
     */
    private Space.Tensor velocity;
    /**
     * Tensor used to sum contributions to velocity dyad, and returned by currentValue() method.
     */
    private Space.Tensor velocityTensor;
    
    public MeterTensorVelocity() {
        this(Simulation.instance);
    }
    public MeterTensorVelocity(Simulation sim) {
        super(sim);
        velocity = sim.space().makeTensor();
        velocityTensor = sim.space().makeTensor();
    }
    
    /**
     * Indicates that this meter does not reference the phase boundary.
     * @return false
     */
    public boolean usesPhaseBoundary() {return false;}
    
    /**
     * Indicates that this meter uses iterators.
     * @return true
     */
    public boolean usesPhaseIteratorFactory() {return true;}
    
    /**
     * This meter needs iterators to do its measurements, so this method overrides the no-op method of AbstractMeter 
     * It obtains the necessary iterators from the phase.
     */
	protected void setPhaseIteratorFactory(IteratorFactory factory) {
        ai1 = factory.makeAtomIterator();
	}
    
    /**
     * Returns the dimension of the measured value, here given as energy
     */
    public Dimension getDimension() {return Dimension.ENERGY;}
    
    /**
     * Descriptive label
     *
     * @return "pp/m"
     */
    public String getLabel() {return "pp/m";}
    
    /**
     * Returns the velocity dyad (mass*vv) summed over all atoms, and divided by N
     */
    public Space.Tensor currentValue() {
        ai1.reset();
        velocityTensor.E(0.0);
        while(ai1.hasNext()) {
            Atom a = ai1.next();
            velocity.E(a.momentum(), a.momentum());
            velocity.TE(a.rm());
            velocityTensor.PE(velocity);
        }
        velocityTensor.TE(1/phase.atomCount());
        return velocityTensor;
    }
    
    /**
     * Returns the velocity dyad (mass*vv) for the given atom.
     */
    public Space.Tensor currentValue(Atom atom) {
        velocity.E(atom.momentum(), atom.momentum());
        velocity.TE(atom.rm());
        return velocity;
    }
}