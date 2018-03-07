/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;

import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;

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
    
    public BoxInflateDeformable(Box box, Space space){
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
        tempTens.E(((BoundaryDeformablePeriodic)box.getBoundary()).getBoundaryTensor());
        tempTensInv.E(tempTens);
        tempTensInv.invert();
        
        
        /*
         * convert the location of each molecule from x,y,z coordinates
         * into coordinates based on the edge vectors, scale, 
         * convert back, and scale the molecule
         */
        Vector translationVector = translator.getTranslationVector();
        // substract 1 from each dimension so that multiplying by it yields
        // the amount each coordinate is to be translated *by* (not to).
        scaleVector.PE(-1.0);

        IMoleculeList molecules = box.getMoleculeList();
        for(int i = 0; i<molecules.size(); i++) {
            IMolecule molecule = molecules.get(i);
            translationVector.E(moleculeCenter.position(molecule));
            tempTensInv.transform(translationVector);
            translationVector.TE(scaleVector);
            tempTens.transform(translationVector);
            groupScaler.actionPerformed(molecule);
        }
        
        //Reverse the subtraction to the scaleVector
        scaleVector.PE(1.0);
//      Then scale the boundary
        for (int i=0; i<scaleVector.getD(); i++) {
            dimVector.E(box.getBoundary().getEdgeVector(i));
            dimVector.TE(scaleVector.getX(i));
            ((BoundaryDeformablePeriodic)box.getBoundary()).setEdgeVector(i, dimVector);
        }
    }

    protected final Tensor tempTens;
    protected final Tensor tempTensInv;
}
