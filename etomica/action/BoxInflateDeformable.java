package etomica.action;

import etomica.api.IAtom;
import etomica.api.IBox;
import etomica.api.IVector;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomPair;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.util.Debug;

/**
 * Performs actions that cause volume of a deformable system to expand, with molecule
 * positions scaled to keep them in the same relative positions. Inflation can
 * be isotropically or anisotropically.
 * 
 * @author cribbin
 */
public class BoxInflateDeformable extends BoxInflate{

    public BoxInflateDeformable(Space space){
        super(space);
        tempTens = space.makeTensor();
        tempTensInv = space.makeTensor();
    }
    
    public BoxInflateDeformable(IBox box, Space space){
        this(space);
        setBox(box);
    }
   
    /**
     * Performs isotropic inflation.
     */
    public void actionPerformed() {
        if(box == null) return;
        
        //First scale the locations of the molecules.
        
        //get the edge vectors, and invert the tensor with them in it.
        tempTensInv.E(((BoundaryDeformablePeriodic)box.getBoundary()).boundaryTensor());
        tempTensInv.invert();
        tempTens.E(((BoundaryDeformablePeriodic)box.getBoundary()).boundaryTensor());        
        
        
        /*
         * convert the location of each molecule from x,y,z coordinates
         * into coordinates based on the edge vectors, scale, 
         * convert back, and scale the molecule
         */
        moleculeIterator.reset();
        IVector translationVector = translator.getTranslationVector();
        // substract 1 from each dimension so that multiplying by it yields
        // the amount each coordinate is to be translated *by* (not to).
        scaleVector.PE(-1.0);
                
        for(IAtom molecule = moleculeIterator.nextAtom(); molecule != null;
                molecule = moleculeIterator.nextAtom()){
            translationVector.E(moleculeCenter.position(molecule));
            tempTensInv.transform(translationVector);
            translationVector.TE(scaleVector);
            tempTens.transform(translationVector);
            groupScaler.actionPerformed(molecule);
        }
        
        //Reverse the subtraction to the scaleVector
        scaleVector.PE(1.0);
//      Then scale the boundary
        IVector dimensions = box.getBoundary().getDimensions();
        dimensions.TE(scaleVector);
        box.setDimensions(dimensions);
    }

    protected Tensor tempTens;
    protected Tensor tempTensInv;

}
