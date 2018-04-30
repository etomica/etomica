/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.action.AtomActionTranslateBy;
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
import etomica.action.AtomActionTransformed;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.space3d.RotationTensor3D;
import etomica.space3d.Tensor3D;
import etomica.util.random.IRandom;
import etomica.util.random.RandomNumberGenerator;

import java.io.Serializable;

/**
 * CoordinateDefinition implementation for nitrogen molecule. The class takes the first
 * space.D values of u to be real space displacements of the molecule center of
 * mass from its nominal position and 2 rotational displacements. 
 * 
 * @author Tai Boon Tan
 */
public class CoordinateDefinitionNitrogen extends CoordinateDefinitionMolecule
        implements Serializable {

	 public CoordinateDefinitionNitrogen(Simulation sim, Box box, Primitive primitive, Basis basis, Space _space) {
		 this(sim, box, primitive, basis, _space, 2);
		
	 }
	
    public CoordinateDefinitionNitrogen(Simulation sim, Box box, Primitive primitive, Basis basis, Space _space, int rotDim) {
    	
    	super(sim, box, primitive, rotDim, basis, _space);
       /*
        * having rotDim is to have separate coordinateDim
        * for the with and without rotational dof of 
        * the beta phase nitrogen
        */
    	this.rotDim = rotDim;
    	rotationTensor = new RotationTensor3D();
    	rotationTensor.E(tensor);
    	
    	xzOrientationTensor = new Tensor[4];
    	yOrientationTensor = new Tensor[4];
    	
    	positionVector = new Vector[4];

    	for(int i=0; i<xzOrientationTensor.length; i++){
    		xzOrientationTensor[i] = space.makeTensor();
    	}
    	
    	for(int i=0; i<yOrientationTensor.length; i++){
    		yOrientationTensor[i] = space.makeTensor();
    	}
    	
    	positionVector[0] = Vector.of(-0.01769046232605003, -0.036932136480007434, -5.9585171577170365E-5);
    	positionVector[1] = Vector.of(0.018633051735877395, 0.037645190581725024, -5.3840682197738796E-5);
    	positionVector[2] = Vector.of(0.01859592800815902, -0.036939857722296486, -6.615619904963137E-5);
    	positionVector[3] = Vector.of(-0.017669323683296188, 0.03766804096641073, -6.834032947314052E-5);
    	
    	axis = space.makeVector();
    	orientVector = space.makeVector();
    	
        orientationManager = new MoleculeAgentManager(sim, box, new OrientationAgentSource());
        atomGroupAction = new MoleculeChildAtomAction(new AtomActionTransformed(lattice.getSpace()));
        
        translateBy = new AtomActionTranslateBy(space);
        atomGroupActionTranslate = new MoleculeChildAtomAction(translateBy); 
        
        random = new RandomNumberGenerator();
    }

    public void initializeCoordinates(int[] nCells) {
        IMoleculeList moleculeList = box.getMoleculeList();

        int basisSize = lattice.getBasis().getScaledCoordinates().length;

        Vector offset = lattice.getSpace().makeVector();
//        Vector[] primitiveVectors = primitive.vectors();
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
        Vector position = lattice.getSpace().makeVector();
        MoleculeArrayList currentList = null;
		
        if (configuration != null){
        	configuration.initializeCoordinates(box);
        }

        //Scale the lattice offset
        if(!isDoLatticeSum){
	        Vector vectorOfMax = space.makeVector();
	        Vector vectorOfMin = space.makeVector();
	        Vector site = space.makeVector();
	        vectorOfMax.E(Double.NEGATIVE_INFINITY);
	        vectorOfMin.E(Double.POSITIVE_INFINITY);
	        
	        while (indexIterator.hasNext()) {
	            site.E((Vector) lattice.site(indexIterator.next()));
	            for (int i=0; i<site.getD(); i++) {
	                vectorOfMax.setX(i, Math.max(site.getX(i),vectorOfMax.getX(i)));
	                vectorOfMin.setX(i, Math.min(site.getX(i),vectorOfMin.getX(i)));
	            }
	        }
	        offset.Ev1Mv2(vectorOfMax, vectorOfMin);
	        offset.TE(-0.5);
	        offset.ME(vectorOfMin);
        }
        
        indexIterator.reset();
        
        /*
         * this initialization only applied to beta phase crystal
         * [isUnitCellA]
         */
        boolean isUnitCellA = true;
        
        for (int iMolecule = 0; iMolecule<moleculeList.size(); iMolecule++) {
            IMolecule molecule = moleculeList.get(iMolecule);
            
            if (isAlpha){
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
            }
            
            if (isGamma){
	            int rotationNum = iMolecule % 2;
	            if (configuration == null) {
	               
	            	// initialize the coordinate
	                molecule.getType().initializeConformation(molecule);
	                
	                // do the orientation
	                ((AtomActionTransformed)atomGroupAction.getAtomAction()).setTransformationTensor(yOrientationTensor[rotationNum]);
	                atomGroupAction.actionPerformed(molecule);	            
	            }
            }
            
            if(isBeta || isBetaHCP){
            	int nCell = (int)Math.pow(moleculeList.size()/2, 1.000001/3.0)/nCells[2];
            	int numMolinZ = 2*nCell;
            	
            	if(iMolecule>0 && iMolecule%numMolinZ == 0){
            		isUnitCellA = !isUnitCellA;
            		
            	}
            	
            	int rotationNum; 
            	if(isUnitCellA){
            		rotationNum = iMolecule % 2;
            		
            		if (configuration == null) {
        	               
            			// initialize the coordinate
       	               	molecule.getType().initializeConformation(molecule);
       	                
       	                // do the orientation
       	                
       	                ((AtomActionTransformed)atomGroupAction.getAtomAction()).setTransformationTensor(yOrientationTensor[rotationNum]);
       	                atomGroupAction.actionPerformed(molecule);
       	                
       	                ((AtomActionTransformed)atomGroupAction.getAtomAction()).setTransformationTensor(xzOrientationTensor[rotationNum]);
       	                atomGroupAction.actionPerformed(molecule);
       	             
       	            }
            	} else {
            		rotationNum = (iMolecule % 2) + 2;
            		
            		if (configuration == null) {
        	               
            			// initialize the coordinate
       	               	molecule.getType().initializeConformation(molecule);
       	                
       	                // do the orientation
       	                ((AtomActionTransformed)atomGroupAction.getAtomAction()).setTransformationTensor(yOrientationTensor[rotationNum]);
       	                atomGroupAction.actionPerformed(molecule);
       	                
       	                ((AtomActionTransformed)atomGroupAction.getAtomAction()).setTransformationTensor(xzOrientationTensor[rotationNum]);
       	                atomGroupAction.actionPerformed(molecule);
       	                
            		}	
            		
            	}
            	
            	
            }
            
            if(isBetaLatticeSum){
            	if(moleculeList.size() != 4){
            		throw new RuntimeException("<CoordinateDefinitionNitrogen> The method has been hard-coded to have 4-molecules ONLY!");
            	}
            		
            	if (configuration == null) {
        	               
            		// initialize the coordinate
       	          	molecule.getType().initializeConformation(molecule);
       	                
       	          	// do the orientation
       	          	((AtomActionTransformed)atomGroupAction.getAtomAction()).setTransformationTensor(yOrientationTensor[iMolecule]);
       	          	atomGroupAction.actionPerformed(molecule);
       	                
       	          	((AtomActionTransformed)atomGroupAction.getAtomAction()).setTransformationTensor(xzOrientationTensor[iMolecule]);
       	          	atomGroupAction.actionPerformed(molecule);
       	             
       	      	}
            	
            }
            
            
            int[] ii = indexIterator.next();
            // ii[0] and ii[1] = unit Cell number
            // ii[2] = molecule number in unit cell
        	//System.out.println(ii[0] +" " + ii[1] + " " + ii[2] + " " + ii[3] );
            
            position.E((Vector)lattice.site(ii));
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
                cells[iCell] = new BasisCell(currentList, lattice.getSpace().makeVector());
                cells[iCell].cellPosition.E(position);
            }
            /*
             * Translate the lattice position of molecules
             * ONLY apply for beta phase nitrogen
             */
            if(isBeta){
            	
            	int rotationNum; 
            	if(isUnitCellA){
            		rotationNum = iMolecule % 2;
            		translateBy.setTranslationVector(positionVector[rotationNum]);
            		atomGroupActionTranslate.actionPerformed(molecule);
	             
            	} else {
            		rotationNum = (iMolecule % 2) + 2;
            		translateBy.setTranslationVector(positionVector[rotationNum]);
            		atomGroupActionTranslate.actionPerformed(molecule);
            	}
            }
            
            if(isBetaLatticeSum){
            	translateBy.setTranslationVector(positionVector[iMolecule]);
        		atomGroupActionTranslate.actionPerformed(molecule);
            }
            
            currentList.add(molecule);
        }
        
        initNominalU(cells[totalCells-1].molecules);

        moleculeSiteManager = new MoleculeAgentManager(sim, box, new MoleculeSiteSource(space, positionDefinition));
        siteManager = new AtomLeafAgentManager<Vector>(new SiteSource(space), box);
    }
    
/*
 * abandoned method!!
 */
	public void setBetaPositionAndOrientation(IMoleculeList molecules){
    	
    	for (int i = 0; i < molecules.size() ; i++){
    		
    		Vector[] orientation = new Vector[3];
    		Vector orientationMol2 = space.makeVector();
    		
    		orientation[0] = space.makeVector();
    		orientation[1] = space.makeVector();
    		orientation[2] = space.makeVector();
    		
    		IMolecule molecule = molecules.get(i);
    		IMolecule molecule2;
    		
    		if(i%2 == 0){
	    		molecule2 = molecules.get(i+1);
    		
    		} else {
    			molecule2 = molecules.get(i-1);
        			
    		}
    		
    	   	Vector molleafPos0 = molecule.getChildList().get(0).getPosition();
    	   	Vector molleafPos1 = molecule.getChildList().get(1).getPosition();
    	 
    	  	Vector mol2leafPos0 = molecule2.getChildList().get(0).getPosition();
    	   	Vector mol2leafPos1 = molecule2.getChildList().get(1).getPosition();
    	   	
    	   	
    	   	
    	   	orientation[0].Ev1Mv2(molleafPos1, molleafPos0);
    	    orientation[0].normalize();
  
    	    orientationMol2.Ev1Mv2(mol2leafPos1, mol2leafPos0);
    	    orientationMol2.normalize();
    	    
    	    if(i%2 == 0){
    	    	orientation[2].E(orientation[0]);
    	    	orientation[2].XE(orientationMol2);
    	    	orientation[2].normalize();
    	    	
    	    } else {
    	    	orientation[2].E(orientationMol2);
    	    	orientation[2].XE(orientation[0]);
    	    	orientation[2].normalize();
    	    	
    	    }
    	    
    	    orientation[1].E(orientation[2]);
    	    orientation[1].XE(orientation[0]);
    	    orientation[1].normalize();
    	    
    	    orientationManager.setAgent(molecule, orientation);
    	    moleculeSiteManager.setAgent(molecule, positionDefinition.position(molecule));
    	   
    	}
    	
    	moleculeSiteManager = new MoleculeAgentManager(sim, box, new MoleculeSiteSource(space, positionDefinition));
        siteManager = new AtomLeafAgentManager<Vector>(new SiteSource(space), box);

    }
    
	public void setNominalReference(IMoleculeList molecules){
    	
    	for (int i = 0; i < molecules.size() ; i++){
    		
    		Vector[] orientation = new Vector[3];
    		
    		orientation[0] = space.makeVector();
    		orientation[1] = space.makeVector();
    		orientation[2] = space.makeVector();
    		
    		IMolecule molecule = molecules.get(i);
    		
    	    orientation[0].E(new double[]{1.0, 0.0, 0.0});
    		orientation[1].E(new double[]{0.0, 1.0, 0.0});
    		orientation[2].E(new double[]{0.0, 0.0, 1.0});
  	    
    	    orientationManager.setAgent(molecule, orientation);
    	    moleculeSiteManager.setAgent(molecule, positionDefinition.position(molecule));
    	   
    	}
    	
    	moleculeSiteManager = new MoleculeAgentManager(sim, box, new MoleculeSiteSource(space, positionDefinition));
        siteManager = new AtomLeafAgentManager<Vector>(new SiteSource(space), box);

    }
	
	
    public void setConfiguration(Configuration configuration){
        this.configuration = configuration;
    }
    
    public void setOrientationVectorGamma(Space space){
    	/*
    	 * Reference : R.L. Mills and A.F. Schuch, PRL 23(20) 1969 pg.1154 Fig1
    	 */
    	rotationTensor.setRotationAxis(Vector.of(0.0, 0.0, 1.0), Math.toRadians(45));
    	yOrientationTensor[0].E(rotationTensor);
    	
    	rotationTensor.setRotationAxis(Vector.of(0.0, 0.0, 1.0), Math.toRadians(-45));
    	yOrientationTensor[1].E(rotationTensor);
  	
    }
    
  
    public void setOrientationVectorAlpha(Space space){
    	/*
    	 * Reference : A. Di Nola et al Acta Cryst. (1970) A26, 144 Fig1
    	 */
    	rotationTensor.setRotationAxis(Vector.of(0.0, 1.0, 0.0), Math.toRadians(-45));
    	yOrientationTensor[0].E(rotationTensor);
    	
    	rotationTensor.setRotationAxis(Vector.of(0.0, 1.0, 0.0), Math.toRadians(-45));
    	yOrientationTensor[1].E(rotationTensor);
    	
    	rotationTensor.setRotationAxis(Vector.of(0.0, 1.0, 0.0), Math.toRadians(45));
    	yOrientationTensor[2].E(rotationTensor);
    	
    	rotationTensor.setRotationAxis(Vector.of(0.0, 1.0, 0.0), Math.toRadians(45));
    	yOrientationTensor[3].E(rotationTensor);
    	
    	/*
    	 * rotation Axis about (-x, 0, z)
    	 * ROTATION ANGLE: arctan(1/sqrt(2)) = 35.26438968deg
    	 */
    	double angle = Math.atan(1/Math.sqrt(2));
    	rotationTensor.setRotationAxis(Vector.of(-1.0 / Math.sqrt(2), 0.0, 1.0 / Math.sqrt(2)), -angle);
    	xzOrientationTensor[0].E(rotationTensor);
    	
    	rotationTensor.setRotationAxis(Vector.of(-1.0 / Math.sqrt(2), 0.0, 1.0 / Math.sqrt(2)), angle);
    	xzOrientationTensor[1].E(rotationTensor);
    
    	/*
    	 * rotation Axis about (-x, 0, -z)
    	 */
       	rotationTensor.setRotationAxis(Vector.of(-1.0 / Math.sqrt(2), 0.0, -1.0 / Math.sqrt(2)),  -angle);
    	xzOrientationTensor[2].E(rotationTensor);
    	
    	rotationTensor.setRotationAxis(Vector.of(-1.0 / Math.sqrt(2), 0.0, -1.0 / Math.sqrt(2)),  angle);
    	xzOrientationTensor[3].E(rotationTensor);
    	
    }
    
    public void setOrientationVectorBetaInitial(Space space){
    	/*
    	 * To determine the minimum energy structure.
    	 */
    	rotationTensor.setRotationAxis(Vector.of(0.0, 1.0, 0.0), Math.toRadians(90.0));
    	yOrientationTensor[0].E(rotationTensor);
    	
    	rotationTensor.setRotationAxis(Vector.of(0.0, 1.0, 0.0), Math.toRadians(90.0));
    	yOrientationTensor[1].E(rotationTensor);
    	
    	rotationTensor.setRotationAxis(Vector.of(0.0, 1.0, 0.0), Math.toRadians( 0.0));
    	yOrientationTensor[2].E(rotationTensor);
    	
    	rotationTensor.setRotationAxis(Vector.of(0.0, 1.0, 0.0), Math.toRadians( 0.0));
    	yOrientationTensor[3].E(rotationTensor);
    	
    	/*
    	 * rotation Axis about (x, -y, 0)
    	 */
    	rotationTensor.setRotationAxis(Vector.of(0.6649708673741139, -0.7468692961581154, 0.0), Math.toRadians(0.0));
    	xzOrientationTensor[0].E(rotationTensor);
    	
    	rotationTensor.setRotationAxis(Vector.of(0.6649136890095793, -0.7469202006691695, 0.0), Math.toRadians(0.0));
    	xzOrientationTensor[1].E(rotationTensor);
    
    	/*
    	 * rotation Axis about (-x, -y, 0)
    	 */
       	rotationTensor.setRotationAxis(Vector.of(-0.6647945167590543, -0.7470262716177362, 0.0),  Math.toRadians(0.0));
    	xzOrientationTensor[2].E(rotationTensor);
    	
    	rotationTensor.setRotationAxis(Vector.of(-0.6649586501964135, -0.7468801734742754, 0.0),  Math.toRadians(0.0));
    	xzOrientationTensor[3].E(rotationTensor);
    	
    }
    
    public void setOrientationVectorBeta(Space space){
    	/*
    	 * Through minimized structure.
    	 * 
    	 * this is probably not convenient and a little confusing in term of the notation.
    	 * yOrientation is actually z-Orientation and
    	 * xzOrientation is actually xyOrientation.
    	 * 
    	 * rotation in clockwise (with rotation axis pointing at you/ out of the page) takes in positive angle
    	 */
    	rotationTensor.setRotationAxis(Vector.of(0.0, 0.0, 1.0), Math.toRadians(41.68009045677981));
    	yOrientationTensor[0].E(rotationTensor);
    	
    	rotationTensor.setRotationAxis(Vector.of(0.0, 0.0, 1.0), Math.toRadians(41.675704190880026));
    	yOrientationTensor[1].E(rotationTensor);
    	
    	rotationTensor.setRotationAxis(Vector.of(0.0, 0.0, 1.0), Math.toRadians(-41.6665632114643));
    	yOrientationTensor[2].E(rotationTensor);
    	
    	rotationTensor.setRotationAxis(Vector.of(0.0, 0.0, 1.0), Math.toRadians(-41.679153227700894));
    	yOrientationTensor[3].E(rotationTensor);
    	
    	/*
    	 * rotation Axis about (x, -y, 0)
    	 */
    	rotationTensor.setRotationAxis(Vector.of(0.6649708673741139, -0.7468692961581154, 0.0), Math.toRadians(-33.44359507638775));
    	xzOrientationTensor[0].E(rotationTensor);
    	
    	rotationTensor.setRotationAxis(Vector.of(0.6649136890095793, -0.7469202006691695, 0.0), Math.toRadians(33.45253957035785));
    	xzOrientationTensor[1].E(rotationTensor);
    
    	/*
    	 * rotation Axis about (-x, -y, 0)
    	 */
       	rotationTensor.setRotationAxis(Vector.of(-0.6647945167590543, -0.7470262716177362, 0.0),  Math.toRadians(33.46313235674357));
    	xzOrientationTensor[2].E(rotationTensor);
    	
    	rotationTensor.setRotationAxis(Vector.of(-0.6649586501964135, -0.7468801734742754, 0.0),  Math.toRadians(-33.447623057331526));
    	xzOrientationTensor[3].E(rotationTensor);
    	
    }
    
    public void setOrientationVectorBetaLatticeSum(Space space, double density, double[][] param){
    	/*
    	 * 
    	 *
    	 */
		    	
    	FindBetaN2AngleFromParameter fromParams = new FindBetaN2AngleFromParameter(space, density, param);
    	double[] alpha = fromParams.getAlpha();
    	double[] beta  = fromParams.getBeta();
    	Vector[] rotationAxis = fromParams.getRotationAxis();
    	Vector[] deviationVector = fromParams.getDeviationVector();
    	
    	for(int i=0; i<4; i++){
    		rotationTensor.setRotationAxis(Vector.of(0.0, 0.0, 1.0), Math.toRadians(alpha[i]));
    		yOrientationTensor[i].E(rotationTensor);
    		
        	rotationTensor.setRotationAxis(rotationAxis[i], Math.toRadians(beta[i]));
        	xzOrientationTensor[i].E(rotationTensor);
        	
        	positionVector[i].PE(deviationVector[i]);
    	}
    	
    }
    
    public Tensor[] getXzOrientationTensor() {
		return xzOrientationTensor;
	}

	public Tensor[] getyOrientationTensor() {
		return yOrientationTensor;
	}

    
    /*
     * 
     */
    
    public double[] calcU(IMoleculeList molecules) {
        
    	super.calcU(molecules);
    	
    	if(rotDim==0){
    		return u;
    	}
        int j = 3;
        
        for (int i = 0; i < molecules.size() ; i++){
        	IMolecule molecule = molecules.get(i);
        	Vector[] siteOrientation = (Vector[])orientationManager.getAgent(molecule);
        	
	    	/*
	    	 * Determine the Orientation of Each Molecule
	    	 */
	    	
	    	Vector leafPos0 = molecule.getChildList().get(0).getPosition();
	    	Vector leafPos1 = molecule.getChildList().get(1).getPosition();
	    	
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
	    	
	    	double l = Math.sqrt(axis.Mv1Squared(siteOrientation[0]));
	    	
	    	if(u4==0.0){
	    		
	    		u[j] = l;
	    		u[j+1] = 0;
	    		
	    	} else {
	            double ratio = Math.abs(u3/u4);

	            u[j+1] =  l/Math.sqrt(ratio*ratio+1);
	    		
                u[j] =  ratio*l/Math.sqrt(ratio*ratio+1);
	    	}

	    	if(u4 < 0.0){
                u[j+1] *= -1;
            }
            if (u3 < 0.0){
                u[j] *= -1;
            }

	    	j += coordinateDim/molecules.size();
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
    	Vector orientationMol2 = space.makeVector();
    	
    	for (int i = 0; i < molecules.size() ; i++){
    		
    		Vector[] orientation = new Vector[3];
    			
    		orientation[0] = space.makeVector();
    		orientation[1] = space.makeVector();
    		orientation[2] = space.makeVector();
    		IMolecule molecule = molecules.get(i);
    		
    	   	/*
    	   	 * Determine the Orientation of Each Molecule Within a basis cell
    	   	 */
    	    	
    	   	Vector leafPos0 = molecule.getChildList().get(0).getPosition();
    	   	Vector leafPos1 = molecule.getChildList().get(1).getPosition();
    	 
    	   	orientation[0].Ev1Mv2(leafPos1, leafPos0);
    	    orientation[0].normalize();
    	    /*
    	     * Refer SetOrientationVector method
    	     * 
    	     * Reassigning the orientation[] to the new molecule postion axis
    	     * 
    	     */
    	    
    	    if(isAlpha){
	    	    if (orientation[0].getX(0) > 0 ){
	    	    	orientation[2].E(new double[]{-orientation[0].getX(0), 0 ,orientation[0].getX(2)});
	    	    	orientation[2].normalize();
	    	    	
	    	    } 
	    	    else {
	    	    	orientation[2].E(new double[]{ orientation[0].getX(0), 0 ,orientation[0].getX(2)});
	    	    	orientation[2].normalize();
	    	    	
	    	    }
    	    }
    	    
    	    if(isGamma){
    	    	if (orientation[0].getX(1) > 0){
    	    		orientation[2].E(new double[]{-orientation[0].getX(0), orientation[0].getX(1), 0.0 });
    	    		orientation[2].normalize();
    	    		
    	    	} else {
    	    		orientation[2].E(new double[]{ orientation[0].getX(0), -orientation[0].getX(1), 0.0 });
    	    		orientation[2].normalize();
    	    	}
    	    }
    	    
    	    if(isBeta || isBetaLatticeSum || isBetaHCP){
        		IMolecule molecule2;
        		
        		if(i%2 == 0){
    	    		molecule2 = molecules.get(i+1);
        		
        		} else {
        			molecule2 = molecules.get(i-1);
            			
        		}
        		
        	  	Vector mol2leafPos0 = molecule2.getChildList().get(0).getPosition();
        	   	Vector mol2leafPos1 = molecule2.getChildList().get(1).getPosition();
        	   	      
        	    orientationMol2.Ev1Mv2(mol2leafPos1, mol2leafPos0);
        	    orientationMol2.normalize();
        	    
        	    if(i%2 == 0){
        	    	orientation[2].E(orientation[0]);
        	    	orientation[2].XE(orientationMol2);
        	    	orientation[2].normalize();
        	    	
        	    } else {
        	    	orientation[2].E(orientationMol2);
        	    	orientation[2].XE(orientation[0]);
        	    	orientation[2].normalize();
        	    	
        	    }
    	    }
    	    
    	    orientation[1].E(orientation[2]);
    	    orientation[1].XE(orientation[0]);
    	    orientation[1].normalize();
    	    
    	    orientationManager.setAgent(molecule, orientation);
    	    	
    	}
    	moleculeSiteManager = new MoleculeAgentManager(sim, box, new MoleculeSiteSource(space, positionDefinition));
        siteManager = new AtomLeafAgentManager<Vector>(new SiteSource(space), box);
    }
    
    public void setIsGamma(){
    	isGamma = true;
    }
    
    public void setIsAlpha(){
    	isAlpha = true;
    }

    public void setIsBeta(){
    	isBeta = true;
    }
    
    public void setIsBetaHCP(){
    	isBetaHCP = true;
    }
    
    public void setIsBetaLatticeSum(){
    	isBetaLatticeSum = true;
    }
    
    public void setIsDoLatticeSum() {
		isDoLatticeSum = true;
	}

	public Vector[] getMoleculeOrientation(IMolecule molecule) {
       /*
        * return the initial Orientation of the molecule
        */
        return (Vector[])orientationManager.getAgent(molecule);
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
    	if(rotDim==2){
	        int j=3;
	        
	        for (int i = 0; i < molecules.size() ; i++){
	        	
	        	IMolecule molecule = molecules.get(i);
	        	
	        	Vector[] siteOrientation = (Vector[])orientationManager.getAgent(molecule);

                double u3 = newU[j];
                double u4 = newU[j+1];
                double check = u3*u3 + u4*u4;
                
	        	//Putting the molecule back to its nominal orientation
	        	// when u3 and u4 equal to zero
	        	if(check==0.0){
	        	    ((ConformationNitrogen)((SpeciesN2)molecule.getType()).getConformation()).initializePositions(molecule.getChildList(), siteOrientation[0]);
	        	    j+= coordinateDim/molecules.size();
	        	    continue;
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
			
	        	/*
		         * a.	
		         */
		    	axis.Ea1Tv1(u3, siteOrientation[1]);
		    	axis.PEa1Tv1(u4, siteOrientation[2]);
		    	axis.normalize();
		    	
		    	/*
		    	 * b.
		    	 */
		    	double x = check*0.5;
		    	double costheta = 1. - x; 
                double sintheta = Math.sqrt((2-x)*x);  // sqrt(1 - (1-x)^2) = sqrt(1 - 1 + 2x - x^2) = sqrt(2x - x^2)
		    	
		    	orientVector.Ea1Tv1(costheta, siteOrientation[0]);
		    	orientVector.PEa1Tv1(sintheta, axis);
		    	if (orientVector.isNaN()) throw new RuntimeException();
		    	((ConformationNitrogen)((SpeciesN2)molecule.getType()).getConformation()).initializePositions(molecule.getChildList(), orientVector);
			    
		    	j += coordinateDim/molecules.size();
		    	
	        }
    	}
        super.setToU(molecules, newU);
    }

    
    public void setToUMoleculei(int moleculei, double[] newU) {
    	
    	if(newU.length != 5){
    		throw new RuntimeException("<CoordinateDefinitionNitrogen> setToUMoleculei method, newU[] length should be 5!");
    	}
    	
		IMolecule molecule = box.getMoleculeList().get(moleculei);
        
        Vector[] siteOrientation = (Vector[])orientationManager.getAgent(molecule);


        double u3 = newU[3];
        double u4 = newU[4];
        double check = u3*u3 + u4*u4;
        
        //Putting the molecule back to its nominal orientation
        // when u3 and u4 equal to zero
        if(check==0.0){
            ((ConformationNitrogen)((SpeciesN2)molecule.getType()).getConformation()).initializePositions(molecule.getChildList(), siteOrientation[0]);
            Vector site = getLatticePosition(molecule);
                		            
            atomActionTranslateTo.setDestination(site);
            atomActionTranslateTo.actionPerformed(molecule);
            
            for (int k = 0; k < site.getD(); k++) {
            	work1.setX(k, site.getX(k) + newU[k]);
            }
    		            
            atomActionTranslateTo.setDestination(work1);
            atomActionTranslateTo.actionPerformed(molecule);
            
            return;
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
         *      by using the components computed in the equation
         *  b. find the rotation axis by crossing vector 'new orientation vector' 
         *      with siteOrientation[0]
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
    
        /*
         * a.   
         */
        axis.Ea1Tv1(u3, siteOrientation[1]);
        axis.PEa1Tv1(u4, siteOrientation[2]);
        axis.normalize();
        
        /*
         * b.
         */
        double x = check*0.5;
        double costheta = 1. - x;    // 
        double sintheta = Math.sqrt((2-x)*x);  // sqrt(1 - (1-x)^2) = sqrt(1 - 1 + 2x - x^2) = sqrt(2x - x^2)
        
        orientVector.Ea1Tv1(costheta, siteOrientation[0]);
        orientVector.PEa1Tv1(sintheta, axis);
        ((ConformationNitrogen)((SpeciesN2)molecule.getType()).getConformation()).initializePositions(molecule.getChildList(), orientVector);
		        
        Vector site = getLatticePosition(molecule);
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
    protected Vector[] positionVector;
	protected final Tensor3D tensor = new Tensor3D(new double [][]{{1.0, 0.0, 0.0},{0.0, 1.0, 0.0},{0.0, 0.0, 1.0}});
    protected final Vector axis, orientVector;
    protected Configuration configuration;
    protected MoleculeAgentManager orientationManager; 
    protected final MoleculeChildAtomAction atomGroupAction;
    protected AtomActionTranslateBy translateBy;
    protected MoleculeChildAtomAction atomGroupActionTranslate;
    public boolean isAlpha=false;
    public boolean isGamma=false;
    public boolean isBeta=false;
    public boolean isBetaHCP=false;
    public boolean isBetaLatticeSum=false;
    public boolean isDoLatticeSum=false;
    protected IRandom random;
    public int rotDim;

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
