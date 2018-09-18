/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;
import etomica.action.AtomAction;
import etomica.atom.IAtom;
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
    	transformationTensor.transform(atom.getPosition());
    }
       
    public Tensor getTransformationTensor(){
    	return transformationTensor;
    }

    /**
     * @param transformationTensor The translation vector to set.  A local copy
     * is made of the given vector.
     */
    public void setTransformationTensor(Tensor transformationTensor) {
        this.transformationTensor.E(transformationTensor);
    }
}
