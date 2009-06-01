package etomica.association;

import etomica.api.IAtom;
import etomica.api.IAtomPositioned;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IRandom;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.space.ISpace;
import etomica.space.IVectorRandom;

public class BiasVolumeSphere extends BiasVolume {
    
    private double radius;
    private double innerRadius;
    private final IVectorRandom work;
    private final IRandom random;
    private IBoundary boundary;
    
    public BiasVolumeSphere(ISpace space, IRandom random){
        super(space);
        this.random = random;
        work = (IVectorRandom)space.makeVector();
        radius = 1.0;//size of the sphere
        innerRadius = 0.9;
    }
    
    public void setBox(IBox box) {
    	boundary = box.getBoundary();
    }
    
    public void setBiasSphereRadius(double radius) {
        this.radius = radius;
    }
    
    public void setBiasSphereInnerRadius(double innerRadius) {
        this.innerRadius = innerRadius;
    }
    
    public double getBiasSphereRadius() {return radius;}
    
    public double biasVolume() {
        double prod = 4.0/3.0*Math.PI*(radius*radius*radius-innerRadius*innerRadius*innerRadius);
        return prod;
    }
    
    
    // Insert atom1 in to the Bonding region of atom2 
    //Bonding region is a cube /rectangle for 2 D
    public void biasInsert(IAtom atom1, IAtom atom2) {
        work.setRandomInSphere(random);
        work.TE(radius);//multiply by radius
        while (work.squared()<innerRadius*innerRadius){
        	work.setRandomInSphere(random);
        	work.TE(radius);
        }
        ((IAtomPositioned) atom1).getPosition().Ev1Pv2(((IAtomPositioned) atom2).getPosition(), work);
    }

    /**
     *Function to check for bonding
     *calculate the distance between atom1 and atom2
     */
    
    public boolean isAssociated(IAtom atom1, IAtom atom2){
    
        work.E(((IAtomPositioned) atom2).getPosition());
        work.ME(((IAtomPositioned) atom1).getPosition());
        boundary.nearestImage(work);
        double r2 = work.squared();
        //System.out.println ("atom1 = "+atom1+ " atom2 = "+atom2);
        //System.out.println ("r2 = "+r2);
        return r2 > innerRadius*innerRadius && r2 < radius*radius;
    }
}
