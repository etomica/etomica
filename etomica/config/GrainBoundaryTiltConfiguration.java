package etomica.config;

import etomica.atom.IAtomLeaf;
import etomica.atom.IAtomPositioned;
import etomica.atom.IMolecule;
import etomica.box.Box;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.IndexIteratorRectangular;
import etomica.space.IVector;
import etomica.space.RotationTensor;
import etomica.species.Species;

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
    Species [] species;
    double cutoff;
    double angle;
    double dist;
    protected Species fixedSpecies, mobileSpecies;
    
    public GrainBoundaryTiltConfiguration(BravaisLatticeCrystal aLatticeTOP, BravaisLatticeCrystal aLatticeBOTTOM, Species [] aSpecies, double aCutoff){
        super();    
        
        latticeTOP = aLatticeTOP;
        latticeBOTTOM = aLatticeBOTTOM;
        cutoff = aCutoff;
        species = aSpecies;
        
        eulerRotationL2BoxTOP = latticeTOP.getSpace().makeRotationTensor();
        eulerRotationB2LatticeTOP = latticeTOP.getSpace().makeRotationTensor();
        eulerRotationL2BoxBOTTOM = latticeBOTTOM.getSpace().makeRotationTensor();
        eulerRotationB2LatticeBOTTOM = latticeBOTTOM.getSpace().makeRotationTensor();
        
    }
    
    public void setFixedSpecies(Species newFixedSpecies) {
        fixedSpecies = newFixedSpecies;
    }
    
    public void setMobileSpecies(Species newMobileSpecies) {
        mobileSpecies = newMobileSpecies;
    }
    
    public Species getFixedSpecies() {
        return fixedSpecies;
    }
    
    public Species getMobileSpecies() {
        return mobileSpecies;
    }
    
    /**
     * Sets tensor for rotation about the indicated axis (0=x,1=y,2=z) by 
     * the given angle. Calls method from RotationTensor3D.
     */
    public void setRotation(int axis, double aAngle){
        angle = aAngle;
        
        eulerRotationL2BoxTOP.setAxial(axis, angle);
        eulerRotationL2BoxBOTTOM.setAxial(axis, -angle);
        
        eulerRotationB2LatticeTOP.E(eulerRotationL2BoxTOP);
        eulerRotationB2LatticeTOP.inverse();
        
        eulerRotationB2LatticeBOTTOM.E(eulerRotationL2BoxBOTTOM);
        eulerRotationB2LatticeBOTTOM.inverse();
                
    }
    
    public void initializeCoordinates(Box box){
        
        /**
         * FILL ATOMS IN TOP DOMAIN
         */
        
        IVector boxCorner = box.getSpace().makeVector();
        IVector latticeVector = box.getSpace().makeVector();
        
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
        int [] ncellsTOP = new int[box.getSpace().D()+1];
        IVector latticeCenter = box.getSpace().makeVector();
        
        for(int i=0; i<ncellsTOP.length-1; i++){
            ncellsTOP[i] = (int)Math.ceil(latticeVector.x(i)*2/latticeTOP.getPrimitive().vectors()[i].x(i));
            latticeCenter.setX(i, 0.5*ncellsTOP[i]*latticeTOP.getPrimitive().vectors()[i].x(i));
        }
        ncellsTOP[ncellsTOP.length-1] = latticeTOP.getBasis().getScaledCoordinates().length;
        
        IndexIteratorRectangular indexIteratorTOP = new IndexIteratorRectangular(box.getSpace().D()+1);
        
        indexIteratorTOP.setSize(ncellsTOP);
        
        
        indexIteratorTOP.reset();
        
        while(indexIteratorTOP.hasNext()){
            
            IVector transformedPosition = box.getSpace().makeVector();
            
            transformedPosition.E((IVector)latticeTOP.site(indexIteratorTOP.next()));
            
            transformedPosition.ME(latticeCenter);
                       
            eulerRotationB2LatticeTOP.transform(transformedPosition);
            
            transformedPosition.setX(2,transformedPosition.x(2)+(0.25*box.getBoundary().getDimensions().x(2)));
            
            // If the atom position is outside the original simulation domain A (top half of simulation box)
            if(!box.getBoundary().getShape().contains(transformedPosition)||transformedPosition.x(2)<0.0){
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
        int [] ncellsBOTTOM = new int[box.getSpace().D()+1];
        
        for(int i=0; i<ncellsBOTTOM.length-1; i++){
            ncellsBOTTOM[i] = (int)Math.ceil(latticeVector.x(i)*2/latticeBOTTOM.getPrimitive().vectors()[i].x(i));
            latticeCenter.setX(i, 0.5*ncellsBOTTOM[i]*latticeBOTTOM.getPrimitive().vectors()[i].x(i));
        }
        ncellsBOTTOM[ncellsBOTTOM.length-1] = latticeBOTTOM.getBasis().getScaledCoordinates().length;
        
        IndexIteratorRectangular indexIteratorBOTTOM = new IndexIteratorRectangular(box.getSpace().D()+1);
        
        indexIteratorBOTTOM.setSize(ncellsBOTTOM);
        
        
        indexIteratorBOTTOM.reset();
        
        while(indexIteratorBOTTOM.hasNext()){
            
            IVector transformedPosition = box.getSpace().makeVector();
            
            transformedPosition.E((IVector)latticeBOTTOM.site(indexIteratorBOTTOM.next()));
            
            transformedPosition.ME(latticeCenter);
                       
            eulerRotationB2LatticeBOTTOM.transform(transformedPosition);
            
            //Notice negative sign for bottom domain
            transformedPosition.setX(2,transformedPosition.x(2)+(-0.25*box.getBoundary().getDimensions().x(2)));
            
            // If the atom position is outside the original simulation domain B (bottom half of simulation box)
            if(!box.getBoundary().getShape().contains(transformedPosition)||transformedPosition.x(2)>0.0){
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
        dist = 0.0001;
        
        IVector rij = box.getSpace().makeVector();
        
        int removeCount = 0;
        double range = 0.0;
        for(int i=0; i<box.getLeafList().getAtomCount()-1; i++){
            for(int j=i+1; j<box.getLeafList().getAtomCount(); j++){
                
                rij.E(((IAtomPositioned)box.getLeafList().getAtom(i)).getPosition());
                rij.ME(((IAtomPositioned)box.getLeafList().getAtom(j)).getPosition());
                
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
