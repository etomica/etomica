/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.paracetamol;

import etomica.action.MoleculeChildAtomAction;
import etomica.atom.AtomLeafAgentManager;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.lattice.IndexIteratorRectangular;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.Primitive;
import etomica.molecule.*;
import etomica.molecule.MoleculeAgentManager.MoleculeAgentSource;
import etomica.normalmode.CoordinateDefinitionMolecule;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;

import java.io.Serializable;

/**
 * CoordinateDefinition implementation for paracetamol molecule. The class takes the first
 * space.D values of u to be real space displacements of the molecule center of
 * mass from its nominal position and 3 rotational displacements. 
 * 
 * @author Andrew Schultz & Tai Tan
 */
public class CoordinateDefinitionParacetamol extends CoordinateDefinitionMolecule
        implements Serializable {

    public CoordinateDefinitionParacetamol(Simulation sim, Box box, Primitive primitive, Basis basis, Space _space) {
    	super(sim, box, primitive, 3, basis, _space);
       
       	axes = new Vector[3];
        axes [0] = space.makeVector();
        axes [1] = space.makeVector();
        axes [2] = space.makeVector();
        com = space.makeVector();
        temp = space.makeVector();
        proj = space.makeVector();
        proja = space.makeVector();
        projb = space.makeVector();
        axisNorm = space.makeVector();
        axisNormPrime = space.makeVector();
        axis0 = space.makeVector();
        axis0Prime = space.makeVector();
        b = space.makeVector();
        bprime = space.makeVector();
        c = space.makeVector();
        
//        x = space.makeVector();
//        y = space.makeVector();
//        z = space.makeVector();
//        xPrime = space.makeVector();
//        zPrime = space.makeVector();
        
        xNorm = space.makeVector();
        yNorm = space.makeVector();
        zNorm = space.makeVector();
        
        yDoublePrime = space.makeVector();
        zDoublePrime = space.makeVector();
        
        xTriplePrime = space.makeVector();
        yTriplePrime = space.makeVector();
        zTriplePrime = space.makeVector();
        zQuadruplePrime = space.makeVector();
        
        rotationL = lattice.getSpace().makeTensor();
        rotationM = lattice.getSpace().makeTensor();
        rotationN = lattice.getSpace().makeTensor();
        
        t = lattice.getSpace().makeTensor();
        
        orientationManager = new MoleculeAgentManager(sim, box, new OrientationAgentSource());
        atomGroupAction = new MoleculeChildAtomAction(new AtomActionTransformed(lattice.getSpace()));
    }

    public void initializeCoordinates(int[] nCells) {
        IMoleculeList moleculeList = box.getMoleculeList();

        int basisSize = lattice.getBasis().getScaledCoordinates().length;

        Vector offset = lattice.getSpace().makeVector();
        Vector[] primitiveVectors = primitive.vectors();
        for (int i=0; i<primitiveVectors.length; i++) {
            offset.PEa1Tv1(nCells[i],primitiveVectors[i]);
        }
        offset.TE(-0.5);
        
        IndexIteratorRectangular indexIterator = new IndexIteratorRectangular(space.D()+1);
        int[] iteratorDimensions = new int[space.D()+1];
        
        
        System.arraycopy(nCells, 0, iteratorDimensions, 0, nCells.length);
        iteratorDimensions[nCells.length] = basisSize;
        indexIterator.setSize(iteratorDimensions);

        int totalCells = 1;
        for (int i=0; i<nCells.length; i++) {
            totalCells *= nCells[i];
        }
        
        cells = new BasisCell[totalCells];
        int iCell = -1;
        // Place molecules
        indexIterator.reset();
        Vector position = lattice.getSpace().makeVector();
        MoleculeArrayList currentList = null;
		if (configuration != null){
        	configuration.initializeCoordinates(box);
        }

        for (int iMolecule = 0; iMolecule<moleculeList.getMoleculeCount(); iMolecule++) {
            IMolecule molecule = moleculeList.getMolecule(iMolecule);
            if (configuration == null) {
                // initialize coordinates of child atoms
                molecule.getType().initializeConformation(molecule);
            }
            
            int[] ii = indexIterator.next();
            
            if(configuration == null){
	            for (int i=0; i<3; i++){
	            	t.setComponent(i, i, basisOrientation[ii[3]][i]);
	            }
	        	
	            /*
	             * Take out these 2 lines for MCMoveHarmonic
	             */
	        	((AtomActionTransformed)atomGroupAction.getAtomAction()).setTransformationTensor(t);
	            atomGroupAction.actionPerformed(molecule);
            }
            
            position.E((Vector)lattice.site(ii));
            position.PE(offset);
            
            /*
             * Take out these 2 lines for MCMoveHarmonic
             */
            if (configuration == null) {
                atomActionTranslateTo.setDestination(position);
                atomActionTranslateTo.actionPerformed(molecule);
            }

            if (ii[space.D()] == 0) {
                if (iCell > -1) {
                    initNominalU(cells[iCell].molecules);
                }
                // new cell
                iCell++;
                currentList = new MoleculeArrayList(basisSize);
                cells[iCell] = new BasisCell(new MoleculeListWrapper(currentList), lattice.getSpace().makeVector());
                cells[iCell].cellPosition.E(position);
            }
            currentList.add(molecule);
        }
        
        initNominalU(cells[totalCells-1].molecules);

        moleculeSiteManager = new MoleculeAgentManager(sim, box, new MoleculeSiteSource(space, positionDefinition));
        siteManager = new AtomLeafAgentManager<Vector>(new SiteSource(space), box);
    }
    
    public void setConfiguration(Configuration configuration){
        this.configuration = configuration;
   }
    
    public void setBasisOrthorhombic(){
    	basisOrientation = new double [8][3];
    	
    	basisOrientation[0][0] =  1;
    	basisOrientation[0][1] =  1;
    	basisOrientation[0][2] =  1;
    	
    	basisOrientation[1][0] = -1;
    	basisOrientation[1][1] = -1;
    	basisOrientation[1][2] =  1;
    	
    	basisOrientation[2][0] = -1;
    	basisOrientation[2][1] =  1;
    	basisOrientation[2][2] = -1;
    	
    	basisOrientation[3][0] =  1;
    	basisOrientation[3][1] = -1;
    	basisOrientation[3][2] = -1;
    	
    	basisOrientation[4][0] = -1;
    	basisOrientation[4][1] = -1;
    	basisOrientation[4][2] = -1;
    	
    	basisOrientation[5][0] =  1;
    	basisOrientation[5][1] =  1;
    	basisOrientation[5][2] = -1;
    	
    	basisOrientation[6][0] =  1;
    	basisOrientation[6][1] = -1;
    	basisOrientation[6][2] =  1;
    	
    	basisOrientation[7][0] = -1;
    	basisOrientation[7][1] =  1;
    	basisOrientation[7][2] =  1;
    	
    }
    
    public void setBasisMonoclinic(){
    	basisOrientation = new double [4][3];
    	
    	basisOrientation[0][0] =  1;
    	basisOrientation[0][1] =  1;
    	basisOrientation[0][2] =  1;
    	
    	basisOrientation[1][0] = -1;
    	basisOrientation[1][1] =  1;
    	basisOrientation[1][2] = -1;
    	
    	basisOrientation[2][0] = -1;
    	basisOrientation[2][1] = -1;
    	basisOrientation[2][2] = -1;
    	
    	basisOrientation[3][0] =  1;
    	basisOrientation[3][1] = -1;
    	basisOrientation[3][2] =  1;
    	
    }
    
    /*
     * 
     */
    
    public double[] calcU(IMoleculeList molecules) {
        
    	super.calcU(molecules);
        int j = 3;
        
        for (int i=0; i < molecules.getMoleculeCount() ; i++){
        	IMolecule molecule = molecules.getMolecule(i);
        	Vector[] siteOrientation = (Vector[])orientationManager.getAgent(molecule);
        	
	    	/*
	    	 * Determine the Orientation of Each Molecule
	    	 */
	    	
	    	Vector leafPos0 = molecule.getChildList().get(0).getPosition();
	    	Vector leafPos5 = molecule.getChildList().get(5).getPosition();
	    	Vector leafPos10 = molecule.getChildList().get(10).getPosition();
	    	
	    	/*
	    	 * Determine axis 1 by using Vector Projection
	    	 */
	    	axis0Prime.Ev1Mv2(leafPos5, leafPos0);
	    	
	    	axisNormPrime.Ev1Mv2(leafPos10, leafPos0);
	    	double dotProd = axisNormPrime.dot(axis0Prime);
	    	proj.Ea1Tv1(dotProd/axis0Prime.squared(),axis0Prime);
	    	axisNormPrime.ME(proj);
	    	
	    	axis0Prime.normalize();
	      	axisNormPrime.normalize();
	    	
	    	u[j] = axis0Prime.dot(siteOrientation[1]);  
	    	u[j+1] = axis0Prime.dot(siteOrientation[2]); 
	      	 
	    	/*
	    	 * Getting rotation angle at x-axis
	    	 * 
	    	 * axisNorm is perpendicular to siteOrientation[0]
	    	 */
        	double cosAlpha = axis0Prime.dot(siteOrientation[0]); 
	    	
        	if (cosAlpha > 0.999999999){
	    		axisNorm.E(axisNormPrime); 
	    		
	    	} else {
	    		
	    		bprime.Ea1Tv1( 1.0 /cosAlpha, siteOrientation[0]);
	    		bprime.ME(axis0Prime);
	    		bprime.normalize();
	    		
	    		b.Ea1Tv1(cosAlpha, bprime);
	    		b.PEa1Tv1(-Math.sqrt(1-cosAlpha*cosAlpha), axis0Prime);
	    		b.normalize();
	    		
	    		double bComponent = axisNormPrime.dot(bprime);
	    		
	    		c.E(axis0Prime);
	    		c.XE(bprime);
	    		double cComponent = axisNormPrime.dot(c);
	    		
	    		axisNorm.Ea1Tv1(bComponent, b);
	    		axisNorm.PEa1Tv1(cComponent, c);
	    		axisNorm.normalize();
	    	}
	    	
	    	double dot = axisNorm.dot(siteOrientation[1]);
	    	
	    	if (dot > 0.999999999999){
	    		u[j+2] = 0;
	    	} else
	    	
	    	if (dot < -0.99999999999){
	    		u[j+2] = Math.PI;
	    		System.out.println("dot<-0.9999");
	    		throw new RuntimeException("less than -1");
	    	} else 
	    	
	    	u[j+2] = Math.acos(dot);
	    	axisNorm.XE(siteOrientation[1]);
	    	axisNorm.normalize();
	    	
	    	if (axisNorm.dot(siteOrientation[0]) < 0){
	    		
		    	u[j+2] = -u[j+2];  
	    	}
	        
	    	j += coordinateDim/molecules.getMoleculeCount();
        }
        return u;
     }

    /**
     * Override if nominal U is more than the lattice position of the molecule
     */
    public void initNominalU(IMoleculeList molecules) {
    	
    	for (int i=0; i < molecules.getMoleculeCount() ; i++){
    		
    		Vector[] orientation = new Vector[3];
    		
    		orientation[0] = space.makeVector();
    		orientation[1] = space.makeVector();
    		orientation[2] = space.makeVector();
    		IMolecule molecule = molecules.getMolecule(i);
    		
    	    	/*
    	    	 * Determine the Orientation of Each Molecule
    	    	 */
    	    	
    	    	Vector leafPos0 = molecule.getChildList().get(0).getPosition();
    	    	Vector leafPos5 = molecule.getChildList().get(5).getPosition();
    	    	Vector leafPos10 = molecule.getChildList().get(10).getPosition();
    	    	
    	    	
    	    	/*
    	    	 * Determine axis 1 by using Vector Projection
    	    	 */
    	    	axes[0].Ev1Mv2(leafPos5, leafPos0);
    	    	axes[1].Ev1Mv2(leafPos10, leafPos0);
    	    	
    	    	double dotProd = axes[0].dot(axes[1]);
    	    	proj.Ea1Tv1(dotProd/axes[0].squared(),axes[0]);
    	    	axes[1].ME(proj);
    	    	
    	    	axes[0].normalize();
    	    	axes[1].normalize();
    	    	
    	    	axes[2].E(axes[0]);
    	    	axes[2].XE(axes[1]);
    	    	axes[2].normalize();
    	    	
    	    	orientation[0].E(axes[0]);
    	    	orientation[1].E(axes[1]);
    	    	orientation[2].E(axes[2]);
    	    	
    	    	orientationManager.setAgent(molecule, orientation);
    	    	
    		}
 
    }

    public void setToU(IMoleculeList molecules, double[] newU) {
    	
        int j=3;
        
        for (int i=0; i < molecules.getMoleculeCount() ; i++){
        	
        	IMolecule molecule = molecules.getMolecule(i);
            Vector[] siteOrientation = (Vector[])orientationManager.getAgent(molecule);
	    	
	    	/*
	    	 *   STEP 1
	    	 * 
	    	 * Determine the Orientation of Each Molecule
	    	 */
            
	    	Vector leafPos0 = molecule.getChildList().get(0).getPosition();
	    	Vector leafPos5 = molecule.getChildList().get(5).getPosition();
	    	Vector leafPos10 = molecule.getChildList().get(10).getPosition();
	    	
	    	/*
	    	 * Determine axis 1 by using Vector Projection
	    	 */
	    	
//	    	xPrime.Ev1Mv2(leafPos5, leafPos0);
//	    	xPrime.normalize();
	    	
//	    	x.E(siteOrientation[0]);
	      
	      	/*
	      	 * finding the tensor that brings the arbiVectorPrime to arbiVector
	      	 */ 
//	      	double xDotxPrime = x.dot(xPrime);
//	      	if (xDotxPrime < 0.999999999){
//	      		
//	      		/*
//	             * Getting the axes for vectorPrime (xPrime, y, zPrime) 
//	             * 
//	             * y is the arbitrary rotation axis
//	             */
//	      		
//	      		y.E(xPrime);
//	            y.getE(x);
//		      	y.normalize();
//		      	
//		      	zPrime.E(xPrime);
//		      	zPrime.getE(y);
//		      	zPrime.normalize();
//		      	
//		      	z.E(x);
//		      	z.getE(y);
//		      	z.normalize();
//	      		
//	      		double x1Prime = xPrime.get(0);
//	      		double x2Prime = xPrime.get(1);
//	      		double x3Prime = xPrime.get(2);
//	      	
//		      	double z1Prime = zPrime.get(0);
//		      	double z2Prime = zPrime.get(1);
//		      	double z3Prime = zPrime.get(2);
//		      	
//		      	double x1 = x.get(0);
//		      	double x2 = x.get(1);
//		      	double x3 = x.get(2);
//		      	double y1 = y.get(0);
//		      	double y2 = y.get(1);
//		      	double y3 = y.get(2);
//		      	double z1 = z.get(0);
//		      	double z2 = z.get(1);
//		      	double z3 = z.get(2);
//		      	
//		      	double L11 = x1Prime*x1 + y1*y1 + z1Prime*z1;
//		      	double L12 = x2Prime*x1 + y2*y1 + z2Prime*z1;
//		      	double L13 = x3Prime*x1 + y3*y1 + z3Prime*z1;
//		      	
//		      	double L21 = x1Prime*x2 + y1*y2 + z1Prime*z2;
//		      	double L22 = x2Prime*x2 + y2*y2 + z2Prime*z2;
//		      	double L23 = x3Prime*x2 + y3*y2 + z3Prime*z2;
//		      	
//		      	double L31 = x1Prime*x3 + y1*y3 + z1Prime*z3;
//		      	double L32 = x2Prime*x3 + y2*y3 + z2Prime*z3;
//		      	double L33 = x3Prime*x3 + y3*y3 + z3Prime*z3;
//		      	
//		      	/*
//		      	double L11 = x.dot(xPrime);
//		      	double L12 = x.dot(y);
//		      	double L13 = x.dot(zPrime);
//		      	
//		      	double L21 = y.dot(xPrime);
//		      	double L22 = y.dot(y);
//		      	double L23 = y.dot(zPrime);
//		      	
//		      	double L31 = z.dot(xPrime);
//		      	double L32 = z.dot(y);
//		      	double L33 = z.dot(zPrime);
//		      	*/
//		    	
//		      	rotationL.setComponent(0, 0, L11);
//		    	rotationL.setComponent(0, 1, L12);
//		    	rotationL.setComponent(0, 2, L13);
//		    	rotationL.setComponent(1, 0, L21);
//		    	rotationL.setComponent(1, 1, L22);
//		    	rotationL.setComponent(1, 2, L23);
//		    	rotationL.setComponent(2, 0, L31);
//		    	rotationL.setComponent(2, 1, L32);
//		    	rotationL.setComponent(2, 2, L33);
//		    	
//		    	if(rotationL.isNaN()){
//		    		System.out.println("RotationL tensor is BAD!");
//		    		System.out.println(rotationL);
//		    		throw new RuntimeException();
//		    	}
//		    	((AtomActionTransformed)atomGroupAction.getAction()).setTransformationTensor(rotationL);
//		        atomGroupAction.actionPerformed(molecule);
//	      	}
	        
	      	/*
	      	 *  STEP 2
	      	 * 
	      	 * Treat x as the new vector along the x-axis (vector formed from atom0 to atom10)
	      	 * and yNorm is normal to x-axis
	      	 */
	      	
	      	
	      	for (int a=0; a<3; a++){
	      		t.setComponent(a, a, basisOrientation[i][a]);
	      	}
            molecule.getType().initializeConformation(molecule);
	      	((AtomActionTransformed)atomGroupAction.getAtomAction()).setTransformationTensor(t);
	      	atomGroupAction.actionPerformed(molecule);

	        
	    	xNorm.Ev1Mv2(leafPos5, leafPos0);
	    	yNorm.Ev1Mv2(leafPos10, leafPos0);
	    	
	    	double dotProd1 = yNorm.dot(xNorm);
	    	proj.Ea1Tv1(dotProd1/xNorm.squared(),xNorm);
	    	yNorm.ME(proj);
	    	
	    	xNorm.normalize();
	    	yNorm.normalize();
	    	
	        zNorm.E(xNorm);
	        zNorm.XE(yNorm);
	        zNorm.normalize();
	        
	        //Getting the DoublePrime axes

	    	yDoublePrime.Ea1Tv1(Math.cos(newU[j+2]),siteOrientation[1]);
	    	yDoublePrime.PEa1Tv1(-Math.sin(newU[j+2]), siteOrientation[2]);
	    	yDoublePrime.normalize();
	        
	        zDoublePrime.E(xNorm);
	        zDoublePrime.XE(yDoublePrime);
	        zDoublePrime.normalize();
	        
	        /*
	      	 * finding the tensor that brings the arbiVector to arbiVectorDoublePrime
	      	 */
	    	
	      	double x1Norm = xNorm.getX(0);
	      	double x2Norm = xNorm.getX(1);
	      	double x3Norm = xNorm.getX(2);
	      	
	      	double y1Norm = yNorm.getX(0);
	      	double y2Norm = yNorm.getX(1);
	      	double y3Norm = yNorm.getX(2);
	      	
	      	double z1Norm = zNorm.getX(0);
	      	double z2Norm = zNorm.getX(1);
	      	double z3Norm = zNorm.getX(2);

	      	double y1DoublePrime = yDoublePrime.getX(0);
	      	double y2DoublePrime = yDoublePrime.getX(1);
	      	double y3DoublePrime = yDoublePrime.getX(2);
	      	double z1DoublePrime = zDoublePrime.getX(0);
	      	double z2DoublePrime = zDoublePrime.getX(1);
	      	double z3DoublePrime = zDoublePrime.getX(2);
	      	
	      	
	      	double M11 = x1Norm*x1Norm + y1Norm*y1DoublePrime + z1Norm*z1DoublePrime;
	      	double M12 = x2Norm*x1Norm + y2Norm*y1DoublePrime + z2Norm*z1DoublePrime;
	      	double M13 = x3Norm*x1Norm + y3Norm*y1DoublePrime + z3Norm*z1DoublePrime;
	      	
	      	double M21 = x1Norm*x2Norm + y1Norm*y2DoublePrime + z1Norm*z2DoublePrime;
	      	double M22 = x2Norm*x2Norm + y2Norm*y2DoublePrime + z2Norm*z2DoublePrime;
	      	double M23 = x3Norm*x2Norm + y3Norm*y2DoublePrime + z3Norm*z2DoublePrime;
	      	
	      	double M31 = x1Norm*x3Norm + y1Norm*y3DoublePrime + z1Norm*z3DoublePrime;
	      	double M32 = x2Norm*x3Norm + y2Norm*y3DoublePrime + z2Norm*z3DoublePrime;
	      	double M33 = x3Norm*x3Norm + y3Norm*y3DoublePrime + z3Norm*z3DoublePrime;
	      	
	        /*
	      	double M11 = xNorm.dot(xNorm);
	      	double M12 = xNorm.dot(yNorm);
	      	double M13 = xNorm.dot(zNorm);
	      	
	      	double M21 = yDoublePrime.dot(xNorm);
	      	double M22 = yDoublePrime.dot(yNorm);
	      	double M23 = yDoublePrime.dot(zNorm);
	      	
	      	double M31 = zDoublePrime.dot(xNorm);
	      	double M32 = zDoublePrime.dot(yNorm);
	      	double M33 = zDoublePrime.dot(zNorm);
	        */
	        
	    	rotationM.setComponent(0, 0, M11);
	    	rotationM.setComponent(0, 1, M12);
	    	rotationM.setComponent(0, 2, M13);
	    	rotationM.setComponent(1, 0, M21);
	    	rotationM.setComponent(1, 1, M22);
	    	rotationM.setComponent(1, 2, M23);
	    	rotationM.setComponent(2, 0, M31);
	    	rotationM.setComponent(2, 1, M32);
	    	rotationM.setComponent(2, 2, M33);
	    	
	    	if(rotationM.isNaN()){
	    		System.out.println("RotationM tensor is BAD!");
	    		System.out.println(rotationM);
	    		throw new RuntimeException();
	    	}
	    	((AtomActionTransformed)atomGroupAction.getAtomAction()).setTransformationTensor(rotationM);
	        atomGroupAction.actionPerformed(molecule);
	    	
	    	/*
	    	 *    STEP  3
	    	 * 
	    	 * 
	    	 *  First we find the component for xTriplePrime
	    	 *  x = sqrt(1 - u[j]^2 - u[j+1]^2)
	    	 */

	        if (Math.abs(newU[j])>1e-9 || Math.abs(newU[j+1])>1e-9){
		        	
		    	double x = Math.sqrt(1 - newU[j]*newU[j] - newU[j+1]*newU[j+1]);
		    	
		    	xTriplePrime.Ea1Tv1(x, siteOrientation[0]);
		    	xTriplePrime.PEa1Tv1(newU[j], siteOrientation[1]);
		    	xTriplePrime.PEa1Tv1(newU[j+1], siteOrientation[2]);
		    	
		    	xTriplePrime.normalize();
		    	
		    	yTriplePrime.E(xTriplePrime);
		    	yTriplePrime.XE(siteOrientation[0]);
		    	yTriplePrime.normalize();
		    
		    	zTriplePrime.E(xTriplePrime);
		    	zTriplePrime.XE(yTriplePrime);
		    	zTriplePrime.normalize();
		    	
		    	zQuadruplePrime.E(xNorm);
		    	zQuadruplePrime.XE(yTriplePrime);
		    	zQuadruplePrime.normalize();
		    	
		      	/*
		      	 * finding the tensor that brings the arbiVectorDoublePrime to arbiVectorTriplePrime
		      	 */
		    	
		      	double x1TriplePrime = xTriplePrime.getX(0);
		      	double x2TriplePrime = xTriplePrime.getX(1);
		      	double x3TriplePrime = xTriplePrime.getX(2);
		    	
		      	double y1TriplePrime = yTriplePrime.getX(0);
		      	double y2TriplePrime = yTriplePrime.getX(1);
		      	double y3TriplePrime = yTriplePrime.getX(2);
		      	
		      	double z1TriplePrime = zTriplePrime.getX(0);
		      	double z2TriplePrime = zTriplePrime.getX(1);
		      	double z3TriplePrime = zTriplePrime.getX(2);
		      	
		      	double z1QuadruplePrime = zQuadruplePrime.getX(0);
		      	double z2QuadruplePrime = zQuadruplePrime.getX(1);
		      	double z3QuadruplePrime = zQuadruplePrime.getX(2);
		      	
		      	double N11 = x1Norm*x1TriplePrime + y1TriplePrime*y1TriplePrime + z1QuadruplePrime*z1TriplePrime;
		      	double N12 = x2Norm*x1TriplePrime + y2TriplePrime*y1TriplePrime + z2QuadruplePrime*z1TriplePrime;
		      	double N13 = x3Norm*x1TriplePrime + y3TriplePrime*y1TriplePrime + z3QuadruplePrime*z1TriplePrime;
		      	
		      	double N21 = x1Norm*x2TriplePrime + y1TriplePrime*y2TriplePrime + z1QuadruplePrime*z2TriplePrime;
		      	double N22 = x2Norm*x2TriplePrime + y2TriplePrime*y2TriplePrime + z2QuadruplePrime*z2TriplePrime;
		      	double N23 = x3Norm*x2TriplePrime + y3TriplePrime*y2TriplePrime + z3QuadruplePrime*z2TriplePrime;
		      	
		      	double N31 = x1Norm*x3TriplePrime + y1TriplePrime*y3TriplePrime + z1QuadruplePrime*z3TriplePrime;
		      	double N32 = x2Norm*x3TriplePrime + y2TriplePrime*y3TriplePrime + z2QuadruplePrime*z3TriplePrime;
		      	double N33 = x3Norm*x3TriplePrime + y3TriplePrime*y3TriplePrime + z3QuadruplePrime*z3TriplePrime;
		      	
		    	
	//	      	double N11 = xTriplePrime.dot(xNorm);
	//	      	double N12 = xTriplePrime.dot(yDoublePrime);
	//	      	double N13 = xTriplePrime.dot(zDoublePrime);
	//	      	
	//	      	double N21 = yTriplePrime.dot(xNorm);
	//	      	double N22 = yTriplePrime.dot(yDoublePrime);
	//	      	double N23 = yTriplePrime.dot(zDoublePrime);
	//	      	
	//	      	double N31 = zTriplePrime.dot(xNorm);
	//	      	double N32 = zTriplePrime.dot(yDoublePrime);
	//	      	double N33 = zTriplePrime.dot(zDoublePrime);
		    	
		    	
		    	rotationN.setComponent(0, 0, N11);
		    	rotationN.setComponent(0, 1, N12);
		    	rotationN.setComponent(0, 2, N13);
		    	rotationN.setComponent(1, 0, N21);
		    	rotationN.setComponent(1, 1, N22);
		    	rotationN.setComponent(1, 2, N23);
		    	rotationN.setComponent(2, 0, N31);
		    	rotationN.setComponent(2, 1, N32);
		    	rotationN.setComponent(2, 2, N33);
		    	
		    	if(rotationN.isNaN()){
		    		System.out.println("RotationN tensor is BAD!");
		    		System.out.println(rotationN);
		    		throw new RuntimeException();
		    	}
		    	((AtomActionTransformed)atomGroupAction.getAtomAction()).setTransformationTensor(rotationN);
		        atomGroupAction.actionPerformed(molecule);
	        }
	    	j += coordinateDim/molecules.getMoleculeCount();
	    	
        }
        super.setToU(molecules, newU);
    }

    private static final long serialVersionUID = 1L;
    protected final Vector[] axes;
    protected double [][] basisOrientation ;
    protected final Vector com, temp, axis0, axis0Prime;
    protected final Vector proj, proja, projb;
    protected final Vector axisNorm, axisNormPrime, b, bprime, c;
//    protected final Vector x, y, z, xPrime, zPrime;
    protected final Vector xNorm, yNorm, zNorm;
    protected final Vector yDoublePrime, zDoublePrime;
    protected final Vector xTriplePrime, yTriplePrime, zTriplePrime, zQuadruplePrime;
    protected final Tensor rotationL, rotationM, rotationN;
    protected Configuration configuration;
    
    
    protected MoleculeAgentManager orientationManager; 
    protected final MoleculeChildAtomAction atomGroupAction;
	private Tensor t;

    protected static class OrientationAgentSource implements MoleculeAgentSource, Serializable {
        
        public OrientationAgentSource() {
        }

        public Object makeAgent(IMolecule atom) {
            return null;
        }
        public void releaseAgent(Object agent, IMolecule atom) {
            //nothing to do
        }
        
        private static final long serialVersionUID = 1L;
    }
    
}
