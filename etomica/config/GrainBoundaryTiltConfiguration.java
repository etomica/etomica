package etomica.config;

import etomica.api.IBox;
import etomica.api.IAtomPositioned;
import etomica.api.IMolecule;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.atom.IAtomLeaf;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.IndexIteratorRectangular;
import etomica.space.Boundary;
import etomica.space.RotationTensor;
import etomica.space.Space;
import etomica.space3d.Vector3D;

/**
 * 
 * 
 * @authors ajschultz, msellers
 */

public class GrainBoundaryTiltConfiguration implements Configuration {
    
    RotationTensor eulerRotationL2BoxTOP;
    RotationTensor eulerRotationB2LatticeTOP;
    RotationTensor eulerRotationL2BoxBOTTOM;
    RotationTensor eulerRotationB2LatticeBOTTOM;
    BravaisLatticeCrystal latticeTOP, latticeBOTTOM;
    ISpecies [] species;
    double cutoff;
    double angle;
    double dist;
    protected ISpecies fixedSpecies, mobileSpecies;
    private Space space;
    IVector [] reciprocal, origin, plane;
    IVector normal;
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
        angle = aAngle;
        
        RotationTensor rotT = latticeTOP.getSpace().makeRotationTensor();
        rotT.setAxial(axis, angle);
        rotT.TE(eulerRotationL2BoxTOP);
        eulerRotationL2BoxTOP.E(rotT);
        eulerRotationB2LatticeTOP.E(eulerRotationL2BoxTOP);
        eulerRotationB2LatticeTOP.invert();
                
    }
    
    public void setRotationBOTTOM(int axis, double aAngle){
        angle = aAngle;
        
        RotationTensor rotT = latticeTOP.getSpace().makeRotationTensor();       
        rotT.setAxial(axis, -angle);
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
    	reciprocal = new IVector[space.D()];
    	origin = new IVector[space.D()];
    	for(int i=0; i<space.D(); i++){
	    	reciprocal[i] = space.makeVector();
	    	origin[i] = space.makeVector();
	    	origin[i].setX(i,1);
	   	}
    	    	
    	reciprocal = latticeTOP.getPrimitive().makeReciprocal().vectors();
    	normal = space.makeVector();
    	for(int i=0; i<space.D(); i++){
    		normal.PEa1Tv1(millerPlane[i], reciprocal[i]);
    	}
    	
    	IVector projection = space.makeVector();
    	//rotate Miller plane into X axis, about Z axis (through XY plane).
    	projection.E(normal);
    	//get XY projection of normal
    	projection.setX(2, 0);
    	double theta = Math.acos( projection.dot(origin[0]) / Math.sqrt(projection.squared()) );
    	setRotationTOP(2,theta);
    	setRotationBOTTOM(2,-theta);
    	
    	//rotate Miller plane into Z axis, about Y axis (through XZ plane).
    	projection.E(normal);
    	//get XZ projection of normal
    	projection.setX(1, 0);
    	double phi = Math.acos( projection.dot(origin[2]) / Math.sqrt(projection.squared()) );
    	setRotationTOP(1,phi);
    	setRotationBOTTOM(1,phi);
    }
    
    /**
     * Resizes simulation box to preserve periodic boundary conditions after rotation of lattice.
     * @param box
     */
    public void setBoxSize(IBox box, int[] boxMultiples){
    	int [] m = new int[millerPlane.length];
    	for(int i=0; i<m.length; i++){
    		m[i] = millerPlane[i];
    		if(m[i]==0){
    			m[i]=1;
    		}
    	}
    	//Create vector array of Miller plane XYZ intercepts (least common multiple at X intercept)
    	plane = new IVector[space.D()];
    	for(int i=0; i<plane.length; i++){
    		plane[i] = space.makeVector();
    	}
    	plane[0].setX(0, m[1]*m[2]*latticeTOP.getLatticeConstants()[0]);
    	plane[1].setX(1, m[0]*m[2]*latticeTOP.getLatticeConstants()[1]);
    	plane[2].setX(2, m[0]*m[1]*latticeTOP.getLatticeConstants()[2]);
    	for(int i=0; i<plane.length; i++){
    		if(millerPlane[i]==0){
    			plane[i].setX(i,0);
    		}
    	}
    	
    	//Creates vector in XY plane pointing to edge of Miller plane X and Y intercepts
    	IVector xyprojection = space.makeVector();
    	xyprojection.setX(0, plane[0].x(0));
    	xyprojection.setX(1, plane[1].x(1));
    	
    	//Find distance from xyprojection to Miller plane Z-intercept
    	IVector xaxisperiod = space.makeVector();
    	xaxisperiod.Ev1Mv2(xyprojection, plane[2]);
    	double xaxispbc = Math.sqrt(xaxisperiod.squared());
    	if(xaxispbc==0){
    		xaxispbc=latticeTOP.getLatticeConstants()[0];
    	}
    	
    	//Find distance from xyprojection to Miller plane X-intercept
    	IVector yaxisperiod = space.makeVector();
    	yaxisperiod.Ev1Mv2(xyprojection, plane[0]);
    	double yaxispbc = Math.sqrt(yaxisperiod.squared());
    	if(yaxispbc==0){
    		yaxispbc=latticeTOP.getLatticeConstants()[1];
    	}
    	
    	//If plane has 3 intercepts (XYZ), must double X distance
    	int dub=0;
    	for(int i=0; i<millerPlane.length; i++){
    		if(millerPlane[i]==0){
    			dub=1;
    		}
    	}
    	if(dub==0){
    		xaxispbc = 2.0*xaxispbc;
    	}
    	
    	//Set Box size
    	box.setDimensions(new Vector3D(xaxispbc*boxMultiples[0], yaxispbc*boxMultiples[1], Math.sqrt(normal.squared())*boxMultiples[2]));
    }
    
    public void initializeCoordinates(IBox box){
        
    	if(!(box.getBoundary() instanceof Boundary)) {
    		throw new RuntimeException("Cannot initialize coordinates for a box containing a non etomica.space.Boundary");
    	}
        /**
         * FILL ATOMS IN TOP DOMAIN
         */
        
        IVector boxCorner = space.makeVector();
        IVector latticeVector = space.makeVector();
        
        // Get extremes of rotated simulation domain A (usually top half of Box)
        for(int i=-1; i<3; i++){
            // Map corners of A domain
            boxCorner.E(box.getBoundary().getDimensions());
            boxCorner.TE(0.5);
            boxCorner.setX(2, boxCorner.x(2)*0.5);
            
            if (i!=-1){
                boxCorner.setX(i, -boxCorner.x(i));
            }
            
            //Transform to lattice space coords
            eulerRotationB2LatticeTOP.transform(boxCorner);
            
            for(int j=0; j<3; j++){
                // Loop over all dimensions at box corner
                if (Math.abs(boxCorner.x(j))>latticeVector.x(j)){
                    latticeVector.setX(j, Math.abs(boxCorner.x(j)));
                }
            } 
        }
        
        //Compute number of cells to fit in box
        int [] ncellsTOP = new int[space.D()+1];
        IVector latticeCenter = space.makeVector();
        
        for(int i=0; i<ncellsTOP.length-1; i++){
            ncellsTOP[i] = (int)Math.ceil(latticeVector.x(i)*2/latticeTOP.getPrimitive().vectors()[i].x(i));
            latticeCenter.setX(i, 0.5*ncellsTOP[i]*latticeTOP.getPrimitive().vectors()[i].x(i));
        }
        ncellsTOP[ncellsTOP.length-1] = latticeTOP.getBasis().getScaledCoordinates().length;
        
        IndexIteratorRectangular indexIteratorTOP = new IndexIteratorRectangular(space.D()+1);
        
        indexIteratorTOP.setSize(ncellsTOP);
        
        
        indexIteratorTOP.reset();
        
        while(indexIteratorTOP.hasNext()){
            
            IVector transformedPosition = space.makeVector();
            
            transformedPosition.E((IVector)latticeTOP.site(indexIteratorTOP.next()));
            
            transformedPosition.ME(latticeCenter);
                       
            eulerRotationB2LatticeTOP.transform(transformedPosition);
            
            transformedPosition.setX(2,transformedPosition.x(2)+(0.25*box.getBoundary().getDimensions().x(2)));
            
            // If the atom position is outside the original simulation domain A (top half of simulation box)
            if(!((Boundary)box.getBoundary()).getShape().contains(transformedPosition)||transformedPosition.x(2)<0.0){
               continue;            
            }
            
            // Check to see if this atom needs to be fixed.
            IMolecule a = null;
            if(transformedPosition.x(2)>(box.getBoundary().getDimensions().x(2)/2.0 - cutoff)){
                a = fixedSpecies.makeMolecule();
            }
            else{
                a = mobileSpecies.makeMolecule();
            }
            box.addMolecule(a);
            ((IAtomPositioned)a.getChildList().getAtom(0)).getPosition().E(transformedPosition);
            
        }
        
        /**
         * FILL ATOMS IN BOTTOM DOMAIN
         */
        // Get extremes of rotated simulation domain B (usually bottom half of Box)
        for(int i=-1; i<3; i++){
            // Map corners of B domain
            boxCorner.E(box.getBoundary().getDimensions());
            boxCorner.TE(0.5);
            boxCorner.setX(2, boxCorner.x(2)*0.5);
            
            if (i!=-1){
                boxCorner.setX(i, -boxCorner.x(i));
            }
            
            //Transform to lattice space coords
            eulerRotationB2LatticeBOTTOM.transform(boxCorner);
            
            for(int j=0; j<3; j++){
                // Loop over all dimensions at box corner
                if (Math.abs(boxCorner.x(j))>latticeVector.x(j)){
                    latticeVector.setX(j, Math.abs(boxCorner.x(j)));
                }
            } 
        }
        
        //Compute number of cells to fit in box
        int [] ncellsBOTTOM = new int[space.D()+1];
        
        for(int i=0; i<ncellsBOTTOM.length-1; i++){
            ncellsBOTTOM[i] = (int)Math.ceil(latticeVector.x(i)*2/latticeBOTTOM.getPrimitive().vectors()[i].x(i));
            latticeCenter.setX(i, 0.5*ncellsBOTTOM[i]*latticeBOTTOM.getPrimitive().vectors()[i].x(i));
        }
        ncellsBOTTOM[ncellsBOTTOM.length-1] = latticeBOTTOM.getBasis().getScaledCoordinates().length;
        
        IndexIteratorRectangular indexIteratorBOTTOM = new IndexIteratorRectangular(space.D()+1);
        
        indexIteratorBOTTOM.setSize(ncellsBOTTOM);
        
        
        indexIteratorBOTTOM.reset();
        
        while(indexIteratorBOTTOM.hasNext()){
            
            IVector transformedPosition = space.makeVector();
            
            transformedPosition.E((IVector)latticeBOTTOM.site(indexIteratorBOTTOM.next()));
            
            transformedPosition.ME(latticeCenter);
                       
            eulerRotationB2LatticeBOTTOM.transform(transformedPosition);
            
            //Notice negative sign for bottom domain
            transformedPosition.setX(2,transformedPosition.x(2)+(-0.25*box.getBoundary().getDimensions().x(2)));
            
            // If the atom position is outside the original simulation domain B (bottom half of simulation box)
            if(!((Boundary)box.getBoundary()).getShape().contains(transformedPosition)||transformedPosition.x(2)>0.0){
               continue;            
            }
            
            // Check to see if this atom needs to be fixed. Notice signs/inequalities
            IMolecule a = null;
            if(transformedPosition.x(2)<(-box.getBoundary().getDimensions().x(2)/2.0 + cutoff)){
                a = fixedSpecies.makeMolecule();
            }
            else{
                a = mobileSpecies.makeMolecule();
            }
            box.addMolecule(a);
            ((IAtomPositioned)a.getChildList().getAtom(0)).getPosition().E(transformedPosition);
            
        }
    
     
        /**
         * REMOVE OVERLAPPING ATOMS AT GRAIN BOUNDARY INTERFACE
         */

        dist = 1.5;
        IVector rij = space.makeVector();
        
        int removeCount = 0;
        double range = 0.0;
        for(int i=0; i<box.getLeafList().getAtomCount()-1; i++){
            for(int j=i+1; j<box.getLeafList().getAtomCount(); j++){
                
                rij.E(((IAtomPositioned)box.getLeafList().getAtom(i)).getPosition());
                rij.ME(((IAtomPositioned)box.getLeafList().getAtom(j)).getPosition());
                box.getBoundary().nearestImage(rij);
                range = rij.squared();
                
                if(range<dist){
                    box.removeMolecule(((IAtomLeaf)box.getLeafList().getAtom(j)).getParentGroup());
                    removeCount++;
               }
            }
        }
        
        //Create gap between grains
        /**
        for(int i=0; i<box.getLeafList().getAtomCount()-1; i++){
             
                
                if(Math.abs(((IAtomPositioned)box.getLeafList().getAtom(i)).getPosition().x(2))<3){
                    box.removeMolecule(box.getLeafList().getAtom(i));
                    removeCount++;
             
            }
        }
        */
        System.out.println("Tilt Grain Boundary of "+angle*180/Math.PI+" degrees created.");
        System.out.println(removeCount+" atoms were within "+Math.sqrt(dist)+" of another, and removed.");
    }
}
