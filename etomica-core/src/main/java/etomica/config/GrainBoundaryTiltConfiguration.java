/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.config;

import etomica.api.ISpecies;
import etomica.box.Box;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.IndexIteratorRectangular;
import etomica.molecule.IMolecule;
import etomica.space.Boundary;
import etomica.space.RotationTensor;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Vector3D;

/**
 * @authors ajschultz, msellers
 */

public class GrainBoundaryTiltConfiguration implements Configuration {
    
    RotationTensor eulerRotationL2BoxTOP;
    RotationTensor eulerRotationB2LatticeTOP;
    RotationTensor eulerRotationL2BoxBOTTOM;
    RotationTensor eulerRotationB2LatticeBOTTOM;
    BravaisLatticeCrystal latticeTOP, latticeBOTTOM;
    ISpecies [] species;
    Vector shiftVector;
    double cutoff;
    double xboxshift, yboxshift;
    double phi, theta;
    double dist;
    double spacing;
    protected ISpecies fixedSpecies, mobileSpecies;
    private Space space;
    Vector[] reciprocal;
    Vector[] origin, plane;
    Vector normal;
    int [] millerPlane;
    
    public GrainBoundaryTiltConfiguration(BravaisLatticeCrystal aLatticeTOP,
    		    BravaisLatticeCrystal aLatticeBOTTOM, ISpecies [] aSpecies,
    		    double aCutoff, Space space){
        super();    
        
        latticeTOP = aLatticeTOP;
        latticeBOTTOM = aLatticeBOTTOM;
        cutoff = aCutoff;
        species = aSpecies;
        this.space = space;
        
        eulerRotationL2BoxTOP = latticeTOP.getSpace().makeRotationTensor();
        eulerRotationB2LatticeTOP = latticeTOP.getSpace().makeRotationTensor();
        eulerRotationL2BoxBOTTOM = latticeBOTTOM.getSpace().makeRotationTensor();
        eulerRotationB2LatticeBOTTOM = latticeBOTTOM.getSpace().makeRotationTensor();
        
        resetRotation();
        
    }
    
    public void setFixedSpecies(ISpecies newFixedSpecies) {
        fixedSpecies = newFixedSpecies;
    }
    
    public void setMobileSpecies(ISpecies newMobileSpecies) {
        mobileSpecies = newMobileSpecies;
    }
    
    public ISpecies getFixedSpecies() {
        return fixedSpecies;
    }
    
    public ISpecies getMobileSpecies() {
        return mobileSpecies;
    }
    
    public void resetRotation(){
    	eulerRotationL2BoxTOP.reset();
    	eulerRotationL2BoxBOTTOM.reset();
    	eulerRotationB2LatticeTOP.reset();
    	eulerRotationB2LatticeBOTTOM.reset();
    }
    
    /**
     * Sets tensor for rotation about the indicated axis (0=x,1=y,2=z) by 
     * the given angle. Calls method from RotationTensor3D.
     */
    public void setRotationTOP(int axis, double aAngle){
        RotationTensor rotT = latticeTOP.getSpace().makeRotationTensor();
        rotT.setAxial(axis, aAngle);
        rotT.TE(eulerRotationL2BoxTOP);
        eulerRotationL2BoxTOP.E(rotT);
        eulerRotationB2LatticeTOP.E(eulerRotationL2BoxTOP);
        eulerRotationB2LatticeTOP.invert();
                
    }
    
    public void setRotationBOTTOM(int axis, double aAngle){
        RotationTensor rotT = latticeBOTTOM.getSpace().makeRotationTensor();       
        rotT.setAxial(axis, aAngle);
        rotT.TE(eulerRotationL2BoxBOTTOM);
        eulerRotationL2BoxBOTTOM.E(rotT);
        eulerRotationB2LatticeBOTTOM.E(eulerRotationL2BoxBOTTOM);
        eulerRotationB2LatticeBOTTOM.invert();
                
    }

    /**
     * Allows a user to specify a plane for GB 
     */
    public void setGBplane(int [] m){
    	millerPlane = m;
    	origin = new Vector[space.D()];
    	for(int i=0; i<space.D(); i++){ 
	    	origin[i] = space.makeVector();
	    	origin[i].setX(i,1);
	   	}
    	    	
    	reciprocal = latticeTOP.getPrimitive().makeReciprocal().vectors();
    	normal = space.makeVector();
    	for(int i=0; i<space.D(); i++){
    		normal.PEa1Tv1(millerPlane[i], reciprocal[i]);
    	}
    	spacing = Math.PI*2.0/Math.sqrt(normal.squared());
    	Vector projection = space.makeVector();
    	
    	//rotate Miller plane into Y axis, about Z axis (through XY plane).
    	projection.E(normal);
    	//get XY projection of normal
    	projection.setX(2, 0.0);
    	theta = Math.acos(projection.dot(origin[1]) / Math.sqrt(projection.squared()));
    	setRotationTOP(2,theta);
    	setRotationBOTTOM(2,theta);
    	
    	//rotate Miller plane into Z axis, about X axis (through YZ plane).
    	projection.E(normal);
    	//get normal after rotation
    	eulerRotationL2BoxTOP.transform(projection);
    	
    	phi = Math.acos(projection.dot(origin[2]) / Math.sqrt(projection.squared()));
    	System.out.println(phi*360/Math.PI);
    	setRotationTOP(0,phi);
    	setRotationBOTTOM(0,-phi);
    }
    
    /**
     * Resizes simulation box to preserve periodic boundary conditions after rotation of lattice.
     * @param box
     */
    public void setBoxSize(Box box, int[] boxMultiples){
    	int [] m = new int[millerPlane.length];
    	for(int i=0; i<m.length; i++){
    		m[i] = millerPlane[i];
    		if(m[i]==0){
    			m[i]=1;
    		}
    	}
    	//Create vector array of Miller plane XYZ intercepts
    	plane = new Vector[space.D()];
    	for(int i=0; i<plane.length; i++){
    		plane[i] = space.makeVector();
    	}
    	plane[0].setX(0, m[1]*latticeTOP.getLatticeConstants()[0]);
    	plane[1].setX(1, m[0]*latticeTOP.getLatticeConstants()[1]);
    	plane[2].setX(2, m[2]*latticeTOP.getLatticeConstants()[2]);
    	for(int i=0; i<plane.length; i++){
    		if(millerPlane[i]==0){
    			plane[i].setX(i,0);
    		}
    	}
    	
    	xboxshift = 0;
    	yboxshift =	0;
    	//Find X periodicity - magnitude of Miller plane intersection of X and Y axis.
    	Vector xaxisperiod = space.makeVector();
    	xaxisperiod.Ev1Mv2(plane[0], plane[1]);
    	double xaxispbc = Math.sqrt(xaxisperiod.squared());
    	xboxshift = (xaxispbc/2.0) - Math.sin((Math.PI/2.0) - theta)*plane[0].getX(0);
    	
    	//Find Y periodicity - magnitude of vector formed by Miler plane intersections of Z axis and XY plane.
    	double yaxispbc = 2.0*plane[2].getX(2)/Math.cos((Math.PI/2.0)-phi);
    	//double yaxispbc = Math.sqrt(Math.pow(2.0*plane[2].x(2),2) + Math.pow(plane[0].x(0),2) + Math.pow(plane[1].x(1),2) );
    	
    	//If plane does not intersect X axis
    	if(millerPlane[0]==0){
    		xboxshift = 0;
    		xaxispbc=latticeTOP.getLatticeConstants()[0];
    	}
    	
    	//If plane does not intersect Y axis
        if(millerPlane[1]==0){
        	plane[0].setX(0, m[2]*latticeTOP.getLatticeConstants()[0]);
        	plane[2].setX(2, m[0]*latticeTOP.getLatticeConstants()[2]);
        	for(int i=0; i<plane.length; i++){
        		if(millerPlane[i]==0){
        			plane[i].setX(i,0);
        		}
        	}        	
        	xaxisperiod.Ev1Mv2(plane[0], plane[2]);
        	yaxispbc = Math.sqrt(xaxisperiod.squared());
            yboxshift = (yaxispbc/2.0) - Math.sin((Math.PI/2.0) - phi)*plane[0].getX(0);
        	
            xboxshift = 0;
        	xaxispbc = latticeTOP.getLatticeConstants()[1];
        }
    	
    	//If plane does not intersect Z axis
    	if(millerPlane[2]==0){
    	    yaxispbc=latticeTOP.getLatticeConstants()[2];
    	}
    	
    	System.out.println(xaxispbc+"    "+yaxispbc);
    	
    	//Set Box size
    	box.getBoundary().setBoxSize(new Vector3D(xaxispbc*boxMultiples[0], yaxispbc*boxMultiples[1], 2.0*Math.PI/Math.sqrt(normal.squared())*boxMultiples[2]));
    }
    
    public void initializeCoordinates(Box box){
        
    	if(!(box.getBoundary() instanceof Boundary)) {
    		throw new RuntimeException("Cannot initialize coordinates for a box containing a non etomica.space.Boundary");
    	}
        /**
         * FILL ATOMS IN TOP DOMAIN
         */
        
        Vector boxCorner = space.makeVector();
        Vector latticeVector = space.makeVector();
        
        // Get extremes of rotated simulation domain A (usually top half of Box)
        for(int i=-1; i<3; i++){
            // Map corners of A domain
            boxCorner.E(box.getBoundary().getBoxSize());
            boxCorner.TE(0.5);
            boxCorner.setX(2, boxCorner.getX(2)*0.5);
            
            if (i!=-1){
                boxCorner.setX(i, -boxCorner.getX(i));
            }
            
            //Transform to lattice space coords
            eulerRotationB2LatticeTOP.transform(boxCorner);
            
            for(int j=0; j<3; j++){
                // Loop over all dimensions at box corner
                if (Math.abs(boxCorner.getX(j))>latticeVector.getX(j)){
                    latticeVector.setX(j, Math.abs(boxCorner.getX(j)));
                }
            } 
        }
        
        //Compute number of cells to fit in box
        int [] ncellsTOP = new int[space.D()+1];
        Vector latticeCenter = space.makeVector();
        
        for(int i=0; i<ncellsTOP.length-1; i++){
            ncellsTOP[i] = (int)Math.ceil(latticeVector.getX(i)*2/latticeTOP.getPrimitive().vectors()[i].getX(i));
            latticeCenter.setX(i, 0.5*ncellsTOP[i]*latticeTOP.getPrimitive().vectors()[i].getX(i));
        }
        
        ncellsTOP[ncellsTOP.length-1] = latticeTOP.getBasis().getScaledCoordinates().length;
        
        IndexIteratorRectangular indexIteratorTOP = new IndexIteratorRectangular(space.D()+1);
        
        indexIteratorTOP.setSize(ncellsTOP);
        
        indexIteratorTOP.reset();
        
        while(indexIteratorTOP.hasNext()){
            
            Vector transformedPosition = space.makeVector();
            
            transformedPosition.E((Vector)latticeTOP.site(indexIteratorTOP.next()));
            
            transformedPosition.ME(latticeCenter);
                       
            eulerRotationL2BoxTOP.transform(transformedPosition);
            
            transformedPosition.setX(2,transformedPosition.getX(2)+(0.25*box.getBoundary().getBoxSize().getX(2)));
            
            
            
            // If the atom position is outside the original simulation domain A (top half of simulation box)
            if(!((Boundary)box.getBoundary()).getShape().contains(transformedPosition)||transformedPosition.getX(2)<-0.0001){
               continue;            
            }
            transformedPosition.PE(shiftVector);
            // Check to see if this atom needs to be fixed.
            IMolecule a = null;
            if(transformedPosition.getX(2)>(box.getBoundary().getBoxSize().getX(2)/2.0 - cutoff)){
                a = fixedSpecies.makeMolecule();
            }
            else{
                a = mobileSpecies.makeMolecule();
            }
            box.addMolecule(a);
            a.getChildList().getAtom(0).getPosition().E(transformedPosition);
            
        }
        
        /**
         * FILL ATOMS IN BOTTOM DOMAIN
         */
        // Get extremes of rotated simulation domain B (usually bottom half of Box)
        for(int i=-1; i<3; i++){
            // Map corners of B domain
            boxCorner.E(box.getBoundary().getBoxSize());
            boxCorner.TE(0.5);
            boxCorner.setX(2, boxCorner.getX(2)*0.5);
            
            if (i!=-1){
                boxCorner.setX(i, -boxCorner.getX(i));
            }
            
            //Transform to lattice space coords
            eulerRotationB2LatticeBOTTOM.transform(boxCorner);
            
            for(int j=0; j<3; j++){
                // Loop over all dimensions at box corner
                if (Math.abs(boxCorner.getX(j))>latticeVector.getX(j)){
                    latticeVector.setX(j, Math.abs(boxCorner.getX(j)));
                }
            } 
        }
        
        //Compute number of cells to fit in box
        int [] ncellsBOTTOM = new int[space.D()+1];
        
        for(int i=0; i<ncellsBOTTOM.length-1; i++){
            ncellsBOTTOM[i] = (int)Math.ceil(latticeVector.getX(i)*2/latticeBOTTOM.getPrimitive().vectors()[i].getX(i));
            latticeCenter.setX(i, 0.5*ncellsBOTTOM[i]*latticeBOTTOM.getPrimitive().vectors()[i].getX(i));
        }
        ncellsBOTTOM[ncellsBOTTOM.length-1] = latticeBOTTOM.getBasis().getScaledCoordinates().length;
        
        IndexIteratorRectangular indexIteratorBOTTOM = new IndexIteratorRectangular(space.D()+1);
        
        indexIteratorBOTTOM.setSize(ncellsBOTTOM);
        
        indexIteratorBOTTOM.reset();
        
        while(indexIteratorBOTTOM.hasNext()){
            
            Vector transformedPosition = space.makeVector();
            
            transformedPosition.E((Vector)latticeBOTTOM.site(indexIteratorBOTTOM.next()));
            
            transformedPosition.ME(latticeCenter);
                       
            eulerRotationL2BoxBOTTOM.transform(transformedPosition);
            
            //Notice negative sign for bottom domain
            transformedPosition.setX(2,transformedPosition.getX(2)+(-0.25*box.getBoundary().getBoxSize().getX(2)));
            
            
            
            // If the atom position is outside the original simulation domain B (bottom half of simulation box)
            if(!((Boundary)box.getBoundary()).getShape().contains(transformedPosition)||transformedPosition.getX(2)>0.0001){
               continue;            
            }
            
            // Check to see if this atom needs to be fixed. Notice signs/inequalities
            IMolecule a = null;
            if(transformedPosition.getX(2)<(-box.getBoundary().getBoxSize().getX(2)/2.0 + cutoff)){
                a = fixedSpecies.makeMolecule();
            }
            else{
                a = mobileSpecies.makeMolecule();
            }
            box.addMolecule(a);
            a.getChildList().getAtom(0).getPosition().E(transformedPosition);
            
        }      
      
        
        /**
         * REMOVE OVERLAPPING ATOMS AT GRAIN BOUNDARY INTERFACE
         */
        int removeCount = 0;
        dist = 0.1;
        Vector rij = space.makeVector();
        
        
        double range = 0.0;
        for(int i=0; i<box.getLeafList().getAtomCount()-1; i++){
            for(int j=i+1; j<box.getLeafList().getAtomCount(); j++){
                
                rij.E(box.getLeafList().getAtom(i).getPosition());
                rij.ME(box.getLeafList().getAtom(j).getPosition());
                box.getBoundary().nearestImage(rij);
                range = rij.squared();
                
                if(range<dist){
                    box.removeMolecule(box.getLeafList().getAtom(j).getParentGroup());
                    removeCount++;
               }
            }
        }
        

        
        if(millerPlane[2]==0){
        	theta = 360*theta/Math.PI;
        	if(theta>90){
        	theta = 180-theta;
        	}
        	System.out.println("Tilt Grain Boundary of "+theta+" degrees created.----");
        }else{
        	phi = 360*phi/Math.PI;
        	if(phi>90){
        	//	phi = 180-phi;
        	}
        	System.out.println("Tilt Grain Boundary of "+phi+" degrees created.");
        }
        
        System.out.println(removeCount+" atoms were within "+Math.sqrt(dist)+" of another, and removed.");
    
    }
}
