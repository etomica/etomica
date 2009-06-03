package etomica.association;

import etomica.api.IAtom;
import etomica.api.IAtomPositioned;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IRandom;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.atom.IAtomOriented;
import etomica.space.ISpace;
import etomica.space.IVectorRandom;

public class BiasVolumeSphereOriented extends BiasVolume {
    
    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private double radius;
    private double innerRadius;
    private final IVectorRandom work;
    private final IRandom random;
    private IBoundary boundary;
    private double ec2;
    
    public BiasVolumeSphereOriented(ISpace space, IRandom random){
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
    	if (true) {
    		throw new RuntimeException ("This method has not been implemented.");
    	}
        double prod = 4.0/3.0*Math.PI*(radius*radius*radius-innerRadius*innerRadius*innerRadius);
        return prod;
    }
    
    
    // Insert atom1 in to the Bonding region of atom2 
    //Bonding region is a cube /rectangle for 2 D
    public void biasInsert(IAtom atom1, IAtom atom2) {
    	if (true) {
    		throw new RuntimeException ("This method has not been implemented.");
    	}
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
//        if (atom1.getLeafIndex() == 500 ||atom2.getLeafIndex() == 500 ||atom1.getLeafIndex() == 501 ||atom2.getLeafIndex() == 501){
//	        System.out.println ("atom1 = "+atom1+ " atom2 = "+atom2);
//	        System.out.println ("r2 = "+r2);
//        }
        if (r2 < innerRadius*innerRadius || r2 > radius*radius) {
        	return false;
        }
        IVector e1 = ((IAtomOriented)atom1).getOrientation().getDirection();
        double er1 = e1.dot(work);
//        if (atom1.getLeafIndex() == 500 ||atom2.getLeafIndex() == 500 ||atom1.getLeafIndex() == 501 ||atom2.getLeafIndex() == 501){
//	        System.out.println ("atom1 = "+atom1+ " atom2 = "+atom2);
//	        System.out.println ("er1 = "+er1);
//        }
        if ( er1 < 0.0 || er1*er1 < ec2*r2) {
        	return false;
        }
        IVector e2 = ((IAtomOriented)atom2).getOrientation().getDirection();
        double er2 = e2.dot(work);
//        if (atom1.getLeafIndex() == 500 ||atom2.getLeafIndex() == 500 ||atom1.getLeafIndex() == 501 ||atom2.getLeafIndex() == 501){
//	        System.out.println ("atom1 = "+atom1+ " atom2 = "+atom2);
//	        System.out.println ("er2 = "+er2);
//        }
        return er2 < 0.0 && er2*er2 > ec2*r2;   
    }
public double getTheta() {return Math.acos(ec2);}
    
    /**
     * Accessor method for angle (in radians) describing width of cone.
     */
    public void setTheta(double t) {
        ec2   = Math.cos(t);
        ec2   = ec2*ec2;
    }
}
