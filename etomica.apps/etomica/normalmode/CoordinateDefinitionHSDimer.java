package etomica.normalmode;

import java.io.Serializable;

import etomica.action.MoleculeChildAtomAction;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.ISimulation;
import etomica.api.IVectorMutable;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.MoleculeAgentManager;
import etomica.atom.MoleculeArrayList;
import etomica.atom.MoleculeListWrapper;
import etomica.atom.MoleculeAgentManager.MoleculeAgentSource;
import etomica.config.Configuration;
import etomica.lattice.IndexIteratorRectangular;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.Primitive;
import etomica.models.nitrogen.AtomActionTransformed;
import etomica.space.ISpace;
import etomica.space.Tensor;
import etomica.space3d.RotationTensor3D;
import etomica.space3d.Tensor3D;
import etomica.units.Degree;

/**
 * CoordinateDefinition implementation for HS Dimer. The class takes the first
 * space.D values of u to be real space displacements of the molecule center of
 * mass from its nominal position and 2 rotational displacements. 
 * 
 * @author Tai Boon Tan
 */
public class CoordinateDefinitionHSDimer extends CoordinateDefinitionMolecule
        implements Serializable {
	
    public CoordinateDefinitionHSDimer(ISimulation sim, IBox box, Primitive primitive, Basis basis, ISpace _space) {
    	
    	super(sim, box, primitive, 2, basis, _space);
       /*
        * having rotDim is to have separate coordinateDim
        * for the with and without rotational dof of 
        * the beta phase nitrogen
        */
    	rotationTensor = new RotationTensor3D();
    	rotationTensor.E(tensor);
    	
    	xzOrientationTensor = new Tensor[4];
    	yOrientationTensor = new Tensor[4];
    	
    	for(int i=0; i<xzOrientationTensor.length; i++){
    		xzOrientationTensor[i] = space.makeTensor();
    	}
    	
    	for(int i=0; i<yOrientationTensor.length; i++){
    		yOrientationTensor[i] = space.makeTensor();
    	}
    	
    	axis = space.makeVector();
    	
        orientationManager = new MoleculeAgentManager(sim, box, new OrientationAgentSource());
        atomGroupAction = new MoleculeChildAtomAction(new AtomActionTransformed(lattice.getSpace()));
        
    }

    public void initializeCoordinates(int[] nCells) {
        IMoleculeList moleculeList = box.getMoleculeList();

        int basisSize = lattice.getBasis().getScaledCoordinates().length;

        IVectorMutable offset = lattice.getSpace().makeVector();
//        IVector[] primitiveVectors = primitive.vectors();
//        for (int i=0; i<primitiveVectors.length; i++) {
//            offset.PEa1Tv1(nCells[i],primitiveVectors[i]);
//        }
//        
//        offset.TE(-0.5);
        
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
        IVectorMutable position = lattice.getSpace().makeVector();
        MoleculeArrayList currentList = null;
		
        if (configuration != null){
        	configuration.initializeCoordinates(box);
        }

        //Scale the lattice offset
        IVectorMutable vectorOfMax = space.makeVector();
        IVectorMutable vectorOfMin = space.makeVector();
        IVectorMutable site = space.makeVector();
        vectorOfMax.E(Double.NEGATIVE_INFINITY);
        vectorOfMin.E(Double.POSITIVE_INFINITY);
        
        while (indexIterator.hasNext()) {
            site.E((IVectorMutable) lattice.site(indexIterator.next()));
            for (int i=0; i<site.getD(); i++) {
                vectorOfMax.setX(i, Math.max(site.getX(i),vectorOfMax.getX(i)));
                vectorOfMin.setX(i, Math.min(site.getX(i),vectorOfMin.getX(i)));
            }
        }
        offset.Ev1Mv2(vectorOfMax, vectorOfMin);
        offset.TE(-0.5);
        offset.ME(vectorOfMin);
        
        indexIterator.reset();
                
        for (int iMolecule = 0; iMolecule<moleculeList.getMoleculeCount(); iMolecule++) {
            IMolecule molecule = moleculeList.getMolecule(iMolecule);
            
            int rotationNum = iMolecule % 4;
            if (configuration == null) {
               
            	// initialize the coordinate
                molecule.getType().initializeConformation(molecule);
                
                // do the orientation
                
                ((AtomActionTransformed)atomGroupAction.getAtomAction()).setTransformationTensor(yOrientationTensor[rotationNum]);
                atomGroupAction.actionPerformed(molecule);
                
                ((AtomActionTransformed)atomGroupAction.getAtomAction()).setTransformationTensor(xzOrientationTensor[rotationNum]);
                atomGroupAction.actionPerformed(molecule);
                
            }
                        
            int[] ii = indexIterator.next();
            // ii[0] and ii[1] = unit Cell number
            // ii[2] = molecule number in unit cell
        	//System.out.println(ii[0] +" " + ii[1] + " " + ii[2] + " " + ii[3] );
            
            position.E((IVectorMutable)lattice.site(ii));
            position.PE(offset);
            
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
        siteManager = new AtomLeafAgentManager(new SiteSource(space), box);
    }
    
    public void setConfiguration(Configuration configuration){
        this.configuration = configuration;
    }
      
    public void setOrientationVectorCP2(){
    	/*
    	 *  Vega, Paras and Monson, JCP 96(12), 9060 (1992)
    	 *  
    	 *  arcsin (Pi/2 - 0.6/ sqrt(3) ) 
    	 */
    	
    	
    	double arcsinValue = Math.PI/2 - Math.asin((1.0 / Math.sqrt(3.0))); 
    		
    	rotationTensor.setRotationAxis(space.makeVector(new double[]{0.0, 1.0, 0.0}), arcsinValue);
    	yOrientationTensor[0].E(rotationTensor);
    	
    	rotationTensor.setRotationAxis(space.makeVector(new double[]{0.0, 1.0, 0.0}), arcsinValue);
    	yOrientationTensor[1].E(rotationTensor);
    	
    	rotationTensor.setRotationAxis(space.makeVector(new double[]{0.0, 1.0, 0.0}), arcsinValue);
    	yOrientationTensor[2].E(rotationTensor);
    	
    	rotationTensor.setRotationAxis(space.makeVector(new double[]{0.0, 1.0, 0.0}), arcsinValue);
    	yOrientationTensor[3].E(rotationTensor);
    	
    	/*
    	 * Rotational Axis at Z-Direction
    	 */
    	rotationTensor.setRotationAxis(space.makeVector(new double[]{0.0, 0.0, 1.0}), -Math.PI/2);
    	xzOrientationTensor[0].E(rotationTensor);
    	
    	rotationTensor.setRotationAxis(space.makeVector(new double[]{0.0, 0.0, 1.0}), -Math.PI/2);
    	xzOrientationTensor[1].E(rotationTensor);
    
       	rotationTensor.setRotationAxis(space.makeVector(new double[]{0.0, 0.0, 1.0}), -Math.PI/2);
    	xzOrientationTensor[2].E(rotationTensor);
    	
    	rotationTensor.setRotationAxis(space.makeVector(new double[]{0.0, 0.0, 1.0}), -Math.PI/2);
    	xzOrientationTensor[3].E(rotationTensor);
    }
    
    public void setRotationAngle(double angle){
    	rotationAngle = angle;
    }
    
    public Tensor[] getXzOrientationTensor() {
		return xzOrientationTensor;
	}

	public Tensor[] getyOrientationTensor() {
		return yOrientationTensor;
	}
    
    public double[] calcU(IMoleculeList molecules) {
        
    	super.calcU(molecules);
    	
        int j = 3;
        
        for (int i=0; i < molecules.getMoleculeCount() ; i++){
        	IMolecule molecule = molecules.getMolecule(i);
        	IVectorMutable [] siteOrientation = (IVectorMutable [])orientationManager.getAgent(molecule);
        	
	    	/*
	    	 * Determine the Orientation of Each Molecule
	    	 */
	    	
	    	IVectorMutable leafPos0 = molecule.getChildList().getAtom(0).getPosition();
	    	IVectorMutable leafPos1 = molecule.getChildList().getAtom(1).getPosition();
	    	
	    	/*
	    	 * Determine u3 and u4 by using Vector Projection
	    	 * - taking the dot-product w.r.t. orientation[1]---y' and orientation[2]---z'
	    	 * 
	    	 * with both u[j] and u[j+1] within the constraints:
	    	 *  a. u3^2 + u4^2 = 2*[ 1 - cos(theta) ]
	    	 *  b. u3/u4 = r.orietation[1]/ r.orientaion[2] = ratio
	    	 *  
	    	 *  
	    	 * 
	    	 */
	    	axis.Ev1Mv2(leafPos1, leafPos0);
	       	axis.normalize();
	    	
	       	double u3 = axis.dot(siteOrientation[1]);  
	    	double u4 = axis.dot(siteOrientation[2]); 
	    	double ratio = Math.abs(u3/u4);
	    	
	    	if(Math.abs(u3)< 1e-12 && Math.abs(u4) < 1e-12){
	    		u[j] = 0.0;
	    		u[j+1] = 0.0;
	    	} else {
	    	
		    	double a = axis.dot(siteOrientation[0]);
		    	double theta = Math.acos(a);
		    	/*
		    	 * 
		    	 */
		    	if(Degree.UNIT.fromSim(theta) > 179.99999){
		    		u[j] = Math.sqrt(2);
		    		u[j+1] = Math.sqrt(2);
		    		
		    	} else {
			    	
			    	if(Math.abs(u4) > -1e-10 && Math.abs(u4) < 1e-10){
			    		
			    		u[j] = Math.sqrt(2*(1-Math.cos(theta)));;
			    		if(u3 < 0.0){
			    			u[j] = -u[j];
			    		}
			    		
			    		u[j+1] = u4;
			    	} else {
			    		if(u4 < 0.0){
			    			u[j+1] = -Math.sqrt(2*(1-Math.cos(theta))/(ratio*ratio+1));
			    		} else {
			    			u[j+1] = Math.sqrt(2*(1-Math.cos(theta))/(ratio*ratio+1));
			    		}
			    		
			    		if (u3 < 0.0){
			    			u[j] = -ratio*Math.sqrt(2*(1-Math.cos(theta))/(ratio*ratio+1));
			    		} else {
			    			u[j] = ratio*Math.sqrt(2*(1-Math.cos(theta))/(ratio*ratio+1));
			    		}
			    	}
		    	}
	    	}
	    	j += coordinateDim/molecules.getMoleculeCount();
        }
        return u;
     }

    /**
     * Override if nominal U is more than the lattice position of the molecule
     * Number of molecules equals to basis atoms
     * 
     * initNomial is the method to define the initial orientation of the molecule
     */
    public void initNominalU(IMoleculeList molecules) {
    	
    	for (int i=0; i < molecules.getMoleculeCount() ; i++){
    		
    		IVectorMutable[] orientation = new IVectorMutable[3]; 
    			
    		orientation[0] = space.makeVector();
    		orientation[1] = space.makeVector();
    		orientation[2] = space.makeVector();
    		IMolecule molecule = molecules.getMolecule(i);
    		
    	   	/*
    	   	 * Determine the Orientation of Each Molecule Within a basis cell
    	   	 */
    	    	
    	   	IVectorMutable leafPos0 = molecule.getChildList().getAtom(0).getPosition();
    	   	IVectorMutable leafPos1 = molecule.getChildList().getAtom(1).getPosition();
    	 
    	   	orientation[0].Ev1Mv2(leafPos1, leafPos0);
    	    orientation[0].normalize();
    	    /*
    	     * Refer SetOrientationVector method
    	     * 
    	     * Reassigning the orientation[] to the new molecule postion axis
    	     * 
    	     */
    	    
    	    if (orientation[0].getX(0) > 0 ){
    	    	orientation[2].E(new double[]{-orientation[0].getX(0), 0 ,orientation[0].getX(2)});
    	    	orientation[2].normalize();
    	    	
    	    } 
    	    else {
    	    	orientation[2].E(new double[]{ orientation[0].getX(0), 0 ,orientation[0].getX(2)});
    	    	orientation[2].normalize();
    	    	
    	    }
    	    
    	    orientation[1].E(orientation[2]);
    	    orientation[1].XE(orientation[0]);
    	    orientation[1].normalize();
    	    
    	    orientationManager.setAgent(molecule, orientation);
    	    	
    	}
    	moleculeSiteManager = new MoleculeAgentManager(sim, box, new MoleculeSiteSource(space, positionDefinition));
        siteManager = new AtomLeafAgentManager(new SiteSource(space), box);
    }
    
	public IVectorMutable[] getMoleculeOrientation(IMolecule molecule) {
       /*
        * return the initial Orientation of the molecule
        */
        return (IVectorMutable[])orientationManager.getAgent(molecule);
    }
    
    public void setToU(IMoleculeList molecules, double[] newU) {
    	
    	/*
    	 * use BoxInflate class to
    	 * move the degrees of freedom for volume fluctuation
    	 * in the last 3 components in u[] array
    	 * 
    	 *  x-direction fluctuation : u[coordinateDim-3]
    	 *  y-direction fluctuation : u[coordinateDim-2]
    	 *  z-direction fluctuation : u[coordinateDim-1]
    	 *  
    	 */
    	//for (int i=0; i<rScale.length; i++){
    	//	double currentDimi = box.getBoundary().getBoxSize().getX(0);
    	//	rScale = initVolume.getX(0)/currentDimi; //rescale the fluctuation to the initial volume
    		
    	//}
    	
//    	inflate.setScale(rScale);
//    	inflate.actionPerformed();
    	int j=3;
	        
        for (int i=0; i < molecules.getMoleculeCount() ; i++){
        	
        	IMolecule molecule = molecules.getMolecule(i);
            IVectorMutable[] siteOrientation = (IVectorMutable[])orientationManager.getAgent(molecule);
	    	
            IVectorMutable rotationAxis = space.makeVector();
	    	RotationTensor3D rotation = new RotationTensor3D();
	    	rotation.E(tensor);
            /*
	    	 *   STEP 1
	    	 * 
	    	 * Determine the Orientation of Each Molecule
	    	 * a. To fine the angle between the molecule orientation and
	    	 * 		orientation[0] (its initial position)
	    	 * b. Take the cross product of the 2 vectors to find its rotation axis
	    	 * c. Use RotationTensor3D to rotate the molecule back to its initial position
	    	 */
         
	    	IVectorMutable leafPos0 = molecule.getChildList().getAtom(0).getPosition();
	    	IVectorMutable leafPos1 = molecule.getChildList().getAtom(1).getPosition();
	    	
	        /*
	         * a.
	         */
	    	axis.Ev1Mv2(leafPos1, leafPos0);
	    	axis.normalize();
	    		    	
	    	double angle = Math.acos(axis.dot(siteOrientation[0]));
	    	
	    	/*
	    	 * b.
	    	 */
	    	
	    	if (Math.abs(angle) > 5e-8){ // make sure we DO NOT cross-product vectors with very small angle
	    		rotationAxis.E(axis);
		    	rotationAxis.XE(siteOrientation[0]);
		    	rotationAxis.normalize();
		    	
		    	/*
		    	 * c. rotating clockwise.
		    	 */
		    	rotation.setRotationAxis(rotationAxis, angle);
		    	
		      	if(rotation.isNaN()){
		    		System.out.println("Step 1 Rotation tensor is BAD!");
		    		System.out.println("Rotation Angle is too small, angle: "+ angle);
		    		System.out.println("Rotation is not necessary");
		    		System.out.println(rotation);
		    		throw new RuntimeException();
		    	}
		        ((AtomActionTransformed)atomGroupAction.getAtomAction()).setTransformationTensor(rotation);
	            atomGroupAction.actionPerformed(molecule);
	    	}
                    
	    	/*
	    	 *    STEP  2
	    	 * 
	    	 * 
	    	 *  First we find the component for siteOrientation[1] by the following equation
	    	 *  x = sqrt(1 - u[j]^2 - u[j+1]^2)  ---eq
	    	 *  All the vectors used are normalized
	    	 *  
	    	 *  a. determine the 'new orientation vector' for the molecule
	    	 *  	by using the components computed in the equation
	    	 *  b. find the rotation axis by crossing vector 'new orientation vector' 
	    	 *  	with siteOrientation[0]
	    	 *  c. rotate the molecule to the given position
	    	 *  
	    	 *  the rotation angle is determine through the equation that satisfies the 
	    	 *  equation below:
	    	 *       u3^2 + u4^2 = 2[ 1- cos(theta) ]
	    	 *  at small theta limit, the equation becomes:
	    	 *       u3^2 + u4^2 = theta^2
	    	 *  
	    	 *  
	    	 */

	    	double u3 = newU[j];
        	double u4 = newU[j+1];
        	double check = u3*u3 + u4*u4;
        	
        	/*
    		 * scale u3 and u4 accordingly so that they will satisfy the
    		 *  condition u3^2 + u4^2 < 4.0
    		 *  
    		 *  Free Rotor
    		 */
        	if((Math.abs(u3) > (Math.sqrt(2)+1e-10) || Math.abs(u4) > (Math.sqrt(2)+1e-10)) 
        			&& (check > 4.0)){
        		System.out.println("FREE ROTOR " + u3 + " " + u4);
        		throw new RuntimeException("<CoordinateDefinitionNitrogen> in setToU method");
//	        		double randU3 = random.nextDouble()*Math.sqrt(2);
//	        		double randU4 = random.nextDouble()*Math.sqrt(2);
//	        		
//	        		u3 = randU3;
//	        		u4 = randU4;
//	        		
//	        		if(newU[j] < 0.0){
//	        			newU[j] = -u3;
//	        		} else {
//	        			newU[j] = u3;
//	        		}
//	        		
//	        		if(newU[j+1] < 0.0){
//	        			newU[j+1] = -u4;
//	        		} else {
//	        			newU[j+1] = u4;
//	        		}
        	}
	  	
  
	        if (Math.abs(newU[j])>1e-7 || Math.abs(newU[j+1])>1e-7){
	        	//System.out.println((j+1) +" " + (j+2) +" perform the rotation");
	        	
	        	/*
		         * a.	
		         */
	        	axis.E(0);
		    	axis.Ea1Tv1(newU[j], siteOrientation[1]);
		    	axis.PEa1Tv1(newU[j+1], siteOrientation[2]);
		    	axis.normalize();
		    	
		    	/*
		    	 * b.
		    	 */
		    	angle = Math.acos(1.0000000000000004 - (newU[j]*newU[j] + newU[j+1]*newU[j+1])*0.5);
		    	if(Math.abs(angle) > 1e-7){
			    	rotationAxis.E(0);
			    	rotationAxis.E(axis);
			    	rotationAxis.XE(siteOrientation[0]);
			    	rotationAxis.normalize();
			    	
			    	rotation.setRotationAxis(rotationAxis, -angle);
			    	
			    	if(rotation.isNaN()){
			    		System.out.println("Step 2 Rotation tensor is BAD!");
			    		System.out.println("Rotation Angle is too small, angle: "+ angle);
			    		System.out.println("Rotation is not necessary");
			    		System.out.println(rotation);
			    		throw new RuntimeException();
			    	}
			    	
			    	((AtomActionTransformed)atomGroupAction.getAtomAction()).setTransformationTensor(rotation);
			        atomGroupAction.actionPerformed(molecule);
		    	}
	        }
	    	j += coordinateDim/molecules.getMoleculeCount();
	    	
        }
        super.setToU(molecules, newU);
    }

    
    public void setToUMoleculei(int moleculei, double[] newU) {
    	
    	if(newU.length != 5){
    		throw new RuntimeException("<CoordinateDefinitionNitrogen> setToUMoleculei method, newU[] length should be 5!");
    	}
    	
		IMolecule molecule = box.getMoleculeList().getMolecule(moleculei);
		IVectorMutable[] siteOrientation = (IVectorMutable[])orientationManager.getAgent(molecule);
		    	
		IVectorMutable rotationAxis = space.makeVector();
		RotationTensor3D rotation = new RotationTensor3D();
		rotation.E(tensor);
	            /*
		    	 *   STEP 1
		    	 * 
		    	 * Determine the Orientation of Each Molecule
		    	 * a. To fine the angle between the molecule orientation and
		    	 * 		orientation[0] (its initial position)
		    	 * b. Take the cross product of the 2 vectors to find its rotation axis
		    	 * c. Use RotationTensor3D to rotate the molecule back to its initial position
		    	 */
	         
		IVectorMutable leafPos0 = molecule.getChildList().getAtom(0).getPosition();
		IVectorMutable leafPos1 = molecule.getChildList().getAtom(1).getPosition();
		    	
		/*
         * a.
         */
    	axis.Ev1Mv2(leafPos1, leafPos0);
    	axis.normalize();
		    		    	
    	double angle = Math.acos(axis.dot(siteOrientation[0]));
		    	
    	/*
    	 * b.
    	 */
		    	
    	if (Math.abs(angle) > 5e-8){ // make sure we DO NOT cross-product vectors with very small angle
    		rotationAxis.E(axis);
	    	rotationAxis.XE(siteOrientation[0]);
	    	rotationAxis.normalize();
			    	
	    	/*
	    	 * c. rotating clockwise.
	    	 */
	    	rotation.setRotationAxis(rotationAxis, angle);
			    	
	      	if(rotation.isNaN()){
	    		System.out.println("Step 1 Rotation tensor is BAD!");
	    		System.out.println("Rotation Angle is too small, angle: "+ angle);
	    		System.out.println("Rotation is not necessary");
	    		System.out.println(rotation);
	    		throw new RuntimeException();
	    	}
	        ((AtomActionTransformed)atomGroupAction.getAtomAction()).setTransformationTensor(rotation);
            atomGroupAction.actionPerformed(molecule);
    	}
	                    
	 	/*	  
	 	 *     STEP  2
		 * 
	   	 * 
	   	 *  First we find the component for siteOrientation[1] by the following equation
	   	 *  x = sqrt(1 - u[j]^2 - u[j+1]^2)  ---eq
	   	 *  All the vectors used are normalized
	   	 *  
	   	 *  a. determine the 'new orientation vector' for the molecule
	   	 *  	by using the components computed in the equation
	   	 *  b. find the rotation axis by crossing vector 'new orientation vector' 
	   	 *  	with siteOrientation[0]
	   	 *  c. rotate the molecule to the given position
	   	 *  
	   	 *  the rotation angle is determine through the equation that satisfies the 
	   	 *  equation below:
	   	 *       u3^2 + u4^2 = 2[ 1- cos(theta) ]
	   	 *  at small theta limit, the equation becomes:
	   	 *       u3^2 + u4^2 = theta^2
	   	 *  
	   	 *  
	   	 */
	
    	double u3 = newU[3];
	    double u4 = newU[4];
	    double check = u3*u3 + u4*u4;
	        	
	   	/*
	   	 * scale u3 and u4 accordingly so that they will satisfy the
	     *  condition u3^2 + u4^2 < 4.0
	     *  
	     *  Free Rotor
	     */
	    if((Math.abs(u3) > (Math.sqrt(2)+1e-10) || Math.abs(u4) > (Math.sqrt(2)+1e-10)) 
	    	&& (check > 4.0)){
	    	System.out.println("FREE ROTOR " + u3 + " " + u4);
	    	throw new RuntimeException("<CoordinateDefinitionNitrogen> in setToU method");
	    }
		  	
	  
        if (Math.abs(newU[3])>1e-7 || Math.abs(newU[4])>1e-7){
        	
        	/*
	         * a.	
	         */
        	axis.E(0);
	    	axis.Ea1Tv1(newU[3], siteOrientation[1]);
	    	axis.PEa1Tv1(newU[4], siteOrientation[2]);
	    	axis.normalize();
			    	
	    	/*
	    	 * b.
	    	 */
	    	angle = Math.acos(1.0000000000000004 - (newU[3]*newU[3] + newU[4]*newU[4])*0.5);
	    	if(Math.abs(angle) > 1e-7){
		    	rotationAxis.E(0);
		    	rotationAxis.E(axis);
		    	rotationAxis.XE(siteOrientation[0]);
		    	rotationAxis.normalize();
				    	
		    	rotation.setRotationAxis(rotationAxis, -angle);
				    	
		    	if(rotation.isNaN()){
		    		System.out.println("Step 2 Rotation tensor is BAD!");
		    		System.out.println("Rotation Angle is too small, angle: "+ angle);
		    		System.out.println("Rotation is not necessary");
		    		System.out.println(rotation);
		    		throw new RuntimeException();
		    	}
				    	
		    	((AtomActionTransformed)atomGroupAction.getAtomAction()).setTransformationTensor(rotation);
		        atomGroupAction.actionPerformed(molecule);
	    	}
        }
		        
        IVectorMutable site = getLatticePosition(molecule);
        for (int k = 0; k < site.getD(); k++) {
        	work1.setX(k, site.getX(k) + newU[k]);
        }
		            
        atomActionTranslateTo.setDestination(work1);
        atomActionTranslateTo.actionPerformed(molecule);
    }
    
    private static final long serialVersionUID = 1L;

    protected final RotationTensor3D rotationTensor;
    protected final Tensor[] xzOrientationTensor;
    protected final Tensor[] yOrientationTensor;
	protected final Tensor3D tensor = new Tensor3D(new double [][]{{1.0, 0.0, 0.0},{0.0, 1.0, 0.0},{0.0, 0.0, 1.0}});
    protected final IVectorMutable axis;
    protected Configuration configuration;
    protected MoleculeAgentManager orientationManager; 
    protected final MoleculeChildAtomAction atomGroupAction;
    protected double rotationAngle;

    protected static class OrientationAgentSource implements MoleculeAgentSource, Serializable {
        
        public OrientationAgentSource() {
        }
        public Class getMoleculeAgentClass() {
            return IVectorMutable [].class;
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