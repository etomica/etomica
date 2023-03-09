/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;

import etomica.box.Box;
import etomica.molecule.CenterOfMass;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Class that perform anisotropic angle fluctuation 
 * 
 * @author Tai Boon Tan
 *
 */
public class BoxInflateAnisotropic extends BoxInflate{
    
    public BoxInflateAnisotropic(Space space){
    	super(space);
    	
        cVector = space.makeVector();
        cVectorOld = space.makeVector();
        
    }
    
    public BoxInflateAnisotropic(Box box, Space space){
    	this(space);
    	setBox(box);
    	deltaX = new double[box.getMoleculeList().size()];
    }
    
    /**
     * Performs anisotropic inflation.
     */
    public void actionPerformed() {
        if(box == null) throw new RuntimeException("oops");
        
        Vector translationVector = box.getSpace().makeVector();
        
        cVectorOld.E(box.getBoundary().getEdgeVector(2));
        
        double cx = cVector.getX(0);
        double deltacx = (cx-cVectorOld.getX(0));
        double cz = cVector.getX(2);
        double slope = deltacx/cz;
        
        IMoleculeList molecules = box.getMoleculeList();
        for(int i = 0; i<molecules.size(); i++) {
            IMolecule molecule = molecules.get(i);
            Vector com = CenterOfMass.position(box, molecule);
            
            // delta_x = slope * z
            double h = com.getX(2);
            deltaX[i] = slope*h;
            translationVector.setX(0, deltaX[i]);
            // atoms will end up unwrapped -- out of the box
            transformMolecule(molecule, com, translationVector);
        }

        // set the edgeVectors according to the scaling before passing it to BoundaryDeformablePeriodic
        // only scale the x-, y- and z-axes
        // for the x-component of edgeVector[2], it is not being scale. This is done so as to fluctuate
        // the beta-angle
        
        ((BoundaryDeformablePeriodic)box.getBoundary()).setEdgeVector(2, cVector);

        // now wrap the atoms back in to the box
        wrapMolecules();
    }
    
    
    public void undo(){
    	IMoleculeList molecules = box.getMoleculeList();
    	
    	for(int i = 0; i<molecules.size(); i++) {
    		IMolecule molecule = molecules.get(i);
            Vector com = CenterOfMass.position(box, molecule);
    		translationVector.E(new double[]{-deltaX[i], 0.0, 0.0});
            transformMolecule(molecule, com, translationVector);
    	}
    	
        ((BoundaryDeformablePeriodic)box.getBoundary()).setEdgeVector(2, cVectorOld);

        wrapMolecules();
    }
        
    public void setCVector(Vector cVec){
    	cVector = cVec;
    }
    
    protected Vector cVector, cVectorOld, translationVector;
    protected double[] deltaX;
	private static final long serialVersionUID = 1L;
}
