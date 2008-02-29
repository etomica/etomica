package etomica.paracetamol;

import java.io.Serializable;

import etomica.action.AtomGroupAction;
import etomica.api.IVector;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomSet;
import etomica.atom.AtomTypeMolecule;
import etomica.atom.AtomsetArrayList;
import etomica.atom.IAtom;
import etomica.atom.IAtomPositioned;
import etomica.atom.IMolecule;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.config.Conformation;
import etomica.lattice.IndexIteratorRectangular;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.Primitive;
import etomica.normalmode.CoordinateDefinitionMolecule;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space3d.IVector3D;

/**
 * CoordinateDefinition implementation for paracetamol molecule. The class takes the first
 * space.D values of u to be real space displacements of the molecule center of
 * mass from its nominal position and 3 rotational displacements. 
 * 
 * @author Andrew Schultz & Tai Tan
 */
public class CoordinateDefinitionParacetamol extends CoordinateDefinitionMolecule
        implements Serializable {

    public CoordinateDefinitionParacetamol(Box box, Primitive primitive, Basis basis, Space _space) {
    	super(box, primitive, 3, basis, _space);
       
       	axes = new IVector3D [3];
        axes [0] = (IVector3D)space.makeVector();
        axes [1] = (IVector3D)space.makeVector();
        axes [2] = (IVector3D)space.makeVector();
        com = (IVector3D)space.makeVector();
        temp = (IVector3D)space.makeVector();
        proj = (IVector3D)space.makeVector();
        proja = (IVector3D)space.makeVector();
        projb = (IVector3D)space.makeVector();
        axisNorm = (IVector3D)space.makeVector();
        axisNormPrime = (IVector3D)space.makeVector();
        axis0 = (IVector3D)space.makeVector();
        axis0Prime = (IVector3D)space.makeVector();
        b = (IVector3D)space.makeVector();
        bprime = (IVector3D)space.makeVector();
        c = (IVector3D)space.makeVector();
        
//        x = (IVector3D)space.makeVector();
//        y = (IVector3D)space.makeVector();
//        z = (IVector3D)space.makeVector();
//        xPrime = (IVector3D)space.makeVector();
//        zPrime = (IVector3D)space.makeVector();
        
        xNorm = (IVector3D)space.makeVector();
        yNorm = (IVector3D)space.makeVector();
        zNorm = (IVector3D)space.makeVector();
        
        yDoublePrime = (IVector3D)space.makeVector();
        zDoublePrime = (IVector3D)space.makeVector();
        
        xTriplePrime = (IVector3D)space.makeVector();
        yTriplePrime = (IVector3D)space.makeVector();
        zTriplePrime = (IVector3D)space.makeVector();
        zQuadruplePrime = (IVector3D)space.makeVector();
        
        rotationL = lattice.getSpace().makeTensor();
        rotationM = lattice.getSpace().makeTensor();
        rotationN = lattice.getSpace().makeTensor();
        
        t = lattice.getSpace().makeTensor();
        
        orientationManager = new AtomAgentManager(new OrientationAgentSource(), box);
        atomGroupAction = new AtomGroupAction(new AtomActionTransformed(lattice.getSpace()));
    }

    public void initializeCoordinates(int[] nCells) {
        AtomSet moleculeList = box.getMoleculeList();

        int basisSize = lattice.getBasis().getScaledCoordinates().length;

        IVector offset = lattice.getSpace().makeVector();
        IVector[] primitiveVectors = primitive.vectors();
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
        IVector position = lattice.getSpace().makeVector();
        AtomArrayList currentList = null;
		if (configuration != null){
        	configuration.initializeCoordinates(box);
        }

        for (int iMolecule = 0; iMolecule<moleculeList.getAtomCount(); iMolecule++) {
            IMolecule molecule = (IMolecule)moleculeList.getAtom(iMolecule);
            if (configuration == null) {
                // initialize coordinates of child atoms
                Conformation config = ((AtomTypeMolecule)molecule.getType()).getConformation();
                config.initializePositions(molecule.getChildList());
            }
            
            int[] ii = indexIterator.next();
            
            if(configuration == null){
	            for (int i=0; i<3; i++){
	            	t.setComponent(i, i, basisOrientation[ii[3]][i]);
	            }
	        	
	            /*
	             * Take out these 2 lines for MCMoveHarmonic
	             */
	        	((AtomActionTransformed)atomGroupAction.getAction()).setTransformationTensor(t);
	            atomGroupAction.actionPerformed(molecule);
            }
            
            position.E((IVector)lattice.site(ii));
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
                currentList = new AtomArrayList(basisSize);
                cells[iCell] = new BasisCell(new AtomsetArrayList(currentList), lattice.getSpace().makeVector());
                cells[iCell].cellPosition.E(position);
            }
            currentList.add(molecule);
        }
        
        initNominalU(cells[totalCells-1].molecules);
        
        siteManager = new AtomAgentManager(new SiteSource(space), box);
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
    
    public double[] calcU(AtomSet molecules) {
        
    	super.calcU(molecules);
        int j = 3;
        
        for (int i=0; i < molecules.getAtomCount() ; i++){
        	IMolecule molecule = (IMolecule)molecules.getAtom(i);
        	IVector3D [] siteOrientation = (IVector3D [])orientationManager.getAgent(molecule);
        	
	    	/*
	    	 * Determine the Orientation of Each Molecule
	    	 */
	    	
	    	IVector leafPos0 = ((IAtomPositioned)molecule.getChildList().getAtom(0)).getPosition();
	    	IVector leafPos5 = ((IAtomPositioned)molecule.getChildList().getAtom(5)).getPosition();
	    	IVector leafPos10 = ((IAtomPositioned)molecule.getChildList().getAtom(10)).getPosition();
	    	
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
	        
	    	j += coordinateDim/molecules.getAtomCount();
        }
        return u;
     }

    /**
     * Override if nominal U is more than the lattice position of the molecule
     */
    public void initNominalU(AtomSet molecules) {
    	
    	for (int i=0; i < molecules.getAtomCount() ; i++){
    		
    		IVector3D[] orientation = new IVector3D[3];
    		
    		orientation[0] = (IVector3D)space.makeVector();
    		orientation[1] = (IVector3D)space.makeVector();
    		orientation[2] = (IVector3D)space.makeVector();
    		IMolecule molecule = (IMolecule)molecules.getAtom(i);
    		
    	    	/*
    	    	 * Determine the Orientation of Each Molecule
    	    	 */
    	    	
    	    	IVector leafPos0 = ((IAtomPositioned)molecule.getChildList().getAtom(0)).getPosition();
    	    	IVector leafPos5 = ((IAtomPositioned)molecule.getChildList().getAtom(5)).getPosition();
    	    	IVector leafPos10 = ((IAtomPositioned)molecule.getChildList().getAtom(10)).getPosition();
    	    	
    	    	
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

    public void setToU(AtomSet molecules, double[] newU) {
    	
        int j=3;
        
        for (int i=0; i < molecules.getAtomCount() ; i++){
        	
        	IMolecule molecule = (IMolecule)molecules.getAtom(i);
            IVector3D[] siteOrientation = (IVector3D[])orientationManager.getAgent(molecule);
	    	
	    	/*
	    	 *   STEP 1
	    	 * 
	    	 * Determine the Orientation of Each Molecule
	    	 */
            
	    	IVector leafPos0 = ((IAtomPositioned)molecule.getChildList().getAtom(0)).getPosition();
	    	IVector leafPos5 = ((IAtomPositioned)molecule.getChildList().getAtom(5)).getPosition();
	    	IVector leafPos10 = ((IAtomPositioned)molecule.getChildList().getAtom(10)).getPosition();
	    	
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
//	            y.XE(x);
//		      	y.normalize();
//		      	
//		      	zPrime.E(xPrime);
//		      	zPrime.XE(y);
//		      	zPrime.normalize();
//		      	
//		      	z.E(x);
//		      	z.XE(y);
//		      	z.normalize();
//	      		
//	      		double x1Prime = xPrime.x(0);
//	      		double x2Prime = xPrime.x(1);
//	      		double x3Prime = xPrime.x(2);
//	      	
//		      	double z1Prime = zPrime.x(0);
//		      	double z2Prime = zPrime.x(1);
//		      	double z3Prime = zPrime.x(2);
//		      	
//		      	double x1 = x.x(0);
//		      	double x2 = x.x(1);
//		      	double x3 = x.x(2);
//		      	double y1 = y.x(0);
//		      	double y2 = y.x(1);
//		      	double y3 = y.x(2);
//		      	double z1 = z.x(0);
//		      	double z2 = z.x(1);
//		      	double z3 = z.x(2);
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
            Conformation config = ((AtomTypeMolecule)molecule.getType()).getConformation();
            config.initializePositions(molecule.getChildList());
	      	((AtomActionTransformed)atomGroupAction.getAction()).setTransformationTensor(t);
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
	    	
	      	double x1Norm = xNorm.x(0);
	      	double x2Norm = xNorm.x(1);
	      	double x3Norm = xNorm.x(2);
	      	
	      	double y1Norm = yNorm.x(0);
	      	double y2Norm = yNorm.x(1);
	      	double y3Norm = yNorm.x(2);
	      	
	      	double z1Norm = zNorm.x(0);
	      	double z2Norm = zNorm.x(1);
	      	double z3Norm = zNorm.x(2);

	      	double y1DoublePrime = yDoublePrime.x(0);
	      	double y2DoublePrime = yDoublePrime.x(1);
	      	double y3DoublePrime = yDoublePrime.x(2);
	      	double z1DoublePrime = zDoublePrime.x(0);
	      	double z2DoublePrime = zDoublePrime.x(1);
	      	double z3DoublePrime = zDoublePrime.x(2);
	      	
	      	
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
	    	((AtomActionTransformed)atomGroupAction.getAction()).setTransformationTensor(rotationM);
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
		    	
		      	double x1TriplePrime = xTriplePrime.x(0);
		      	double x2TriplePrime = xTriplePrime.x(1);
		      	double x3TriplePrime = xTriplePrime.x(2);
		    	
		      	double y1TriplePrime = yTriplePrime.x(0);
		      	double y2TriplePrime = yTriplePrime.x(1);
		      	double y3TriplePrime = yTriplePrime.x(2);
		      	
		      	double z1TriplePrime = zTriplePrime.x(0);
		      	double z2TriplePrime = zTriplePrime.x(1);
		      	double z3TriplePrime = zTriplePrime.x(2);
		      	
		      	double z1QuadruplePrime = zQuadruplePrime.x(0);
		      	double z2QuadruplePrime = zQuadruplePrime.x(1);
		      	double z3QuadruplePrime = zQuadruplePrime.x(2);
		      	
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
		    	((AtomActionTransformed)atomGroupAction.getAction()).setTransformationTensor(rotationN);
		        atomGroupAction.actionPerformed(molecule);
	        }
	    	j += coordinateDim/molecules.getAtomCount();
	    	
        }
        super.setToU(molecules, newU);
    }

    private static final long serialVersionUID = 1L;
    protected final IVector3D [] axes;
    protected double [][] basisOrientation ;
    protected final IVector3D com, temp, axis0, axis0Prime; 
    protected final IVector3D proj, proja, projb;
    protected final IVector3D axisNorm, axisNormPrime, b, bprime, c;
//    protected final IVector3D x, y, z, xPrime, zPrime;
    protected final IVector3D xNorm, yNorm, zNorm;
    protected final IVector3D yDoublePrime, zDoublePrime;
    protected final IVector3D xTriplePrime, yTriplePrime, zTriplePrime, zQuadruplePrime;
    protected final Tensor rotationL, rotationM, rotationN;
    protected Configuration configuration;
    
    
    protected AtomAgentManager orientationManager; 
    protected final AtomGroupAction atomGroupAction;
	private Tensor t;

    protected static class OrientationAgentSource implements AgentSource, Serializable {
        
        public OrientationAgentSource() {
        }
        public Class getAgentClass() {
            return IVector [].class;
        }
        public Object makeAgent(IAtom atom) {
            return null;
        }
        public void releaseAgent(Object agent, IAtom atom) {
            //nothing to do
        }
        
        private static final long serialVersionUID = 1L;
    }
    
}