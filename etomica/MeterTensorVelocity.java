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
    private final AtomIteratorList ai1 = new AtomIteratorList();
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
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Velocity tensor, formed from averaging dyad of velocity vector for each atom");
        return info;
    }
    
    /**
     * This meter needs iterators to do its measurements, so this method overrides the no-op method of AbstractMeter 
     * It obtains the necessary iterators from the phase.
     */
	public void setPhase(Phase p) {
	    super.setPhase(p);
        ai1.setBasis(p.speciesMaster.atomList);
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
            velocity.E(a.coord.momentum(), a.coord.momentum());
            velocity.TE(a.coord.rm());
            velocityTensor.PE(velocity);
        }
        velocityTensor.TE(1/phase.atomCount());
        return velocityTensor;
    }
    
    /**
     * Returns the velocity dyad (mass*vv) for the given atom.
     */
    public Space.Tensor currentValue(Atom atom) {
        velocity.E(atom.coord.momentum(), atom.coord.momentum());
        velocity.TE(atom.coord.rm());
        return velocity;
    }
}