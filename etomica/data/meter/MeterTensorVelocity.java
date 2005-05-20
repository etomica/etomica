package etomica.data.meter;
import etomica.Atom;
import etomica.AtomTypeLeaf;
import etomica.EtomicaInfo;
import etomica.Phase;
import etomica.Space;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.atom.iterator.AtomIteratorPhaseDependent;
import etomica.space.ICoordinateKinetic;
import etomica.space.Tensor;
import etomica.units.Dimension;

/**
 * A meter to compute the velocity component of the pressure tensor. 
 * Averages a tensor quantity formed from a dyad of the velocity of each atom. 
 * Specifically, the quantity averaged is 1/N * sum(pp/m), where p is the momentum,
 * m is the mass, and the sum is over all N atoms.
 * 
 * @author Rob Riggleman
 */

public class MeterTensorVelocity extends MeterTensor /*implements MeterTensor.Atomic*/ {
    /**
     * Iterator of atoms.
     */
    private final AtomIteratorPhaseDependent ai1 = new AtomIteratorLeafAtoms();
    /**
     * Tensor used to form velocity dyad for each atom, and returned by currentValue(atom) method.
     */
    private Tensor velocity;
    /**
     * Tensor used to sum contributions to velocity dyad, and returned by currentValue() method.
     */
    private Tensor velocityTensor;
    
    public MeterTensorVelocity(Space space) {
        super(space);
        velocity = space.makeTensor();
        velocityTensor = space.makeTensor();
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Velocity tensor, formed from averaging dyad of velocity vector for each atom");
        return info;
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
    public Tensor getDataAsTensor(Phase p) {
        ai1.setPhase(p);
        ai1.reset();
        velocityTensor.E(0.0);
        int count = 0;
        while(ai1.hasNext()) {
            Atom a = ai1.nextAtom();
            velocity.E(((ICoordinateKinetic)a.coord).velocity(), ((ICoordinateKinetic)a.coord).velocity());
            velocity.TE(((AtomTypeLeaf)a.type).rm());
            velocityTensor.PE(velocity);
            count++;
        }
        velocityTensor.TE(1.0/count);
        return velocityTensor;
    }
    
    /**
     * Returns the velocity dyad (mass*vv) for the given atom.
     */
//    public Space.Tensor currentValue(Atom atom) {
//        velocity.E(atom.coord.momentum(), atom.coord.momentum());
//        velocity.TE(atom.coord.rm());
//        return velocity;
//    }
}
