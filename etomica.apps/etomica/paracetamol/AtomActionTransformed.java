package etomica.paracetamol;
import etomica.action.AtomAction;
import etomica.atom.IAtom;
import etomica.atom.IAtomPositioned;
import etomica.space.Space;
import etomica.space.Tensor;

/**
 * 
 * Moves (translates) an atom by a specified vector amount.
 * To move all atoms in a molecule (or atom group), wrap an
 * instance of this class in an AtomGroupAction.
 * 
 * @author David Kofke
 */
public class AtomActionTransformed implements AtomAction {
    
    private static final long serialVersionUID = 1L;
    private final Tensor transformationTensor;
    
    public AtomActionTransformed(Space space) {
        transformationTensor = space.makeTensor();
    }
    
    public void actionPerformed(IAtom atom) {
    	transformationTensor.transform(((IAtomPositioned)atom).getPosition());
    }
       
    public Tensor getTransformationTensor(){
    	return transformationTensor;
    }

    /**
     * @param destination The translation vector to set.  A local copy
     * is made of the given vector.
     */
    public void setTransformationTensor(Tensor transformationTensor) {
        this.transformationTensor.E(transformationTensor);
    }
}