package simulate;

import java.beans.*;
import java.awt.*;

public final class PotentialHeavyHardDisk extends PotentialHardDisk
{

    public PotentialHeavyHardDisk(double d) {
        super(d);
        setCollisionDiameter(d);
    }

    public double collisionTime(AtomHard atom1, AtomHard atom2) {
        parentPhase.space.uEr1Mr2(r12,atom2.r,atom1.r);  //use instance method   //r2-r1
        Space.uEa1Tv1Ma2Tv2(v12,atom2.rm,atom2.p,atom1.rm,atom1.p);  //v2-v1 = (p/m)2 - (p/m)1
        double bij = Space.v1Dv2(r12,v12);                           //r12 . v12
        double time = Double.MAX_VALUE;

        if(bij < 0.0) {
          double r2 = Space.v1S(r12);
          double velocitySquared = Space.v1S(v12);
          double discriminant = bij*bij - velocitySquared * ( r2 - sig2 );
          if(discriminant > 0) {
            time = (-bij - Math.sqrt(discriminant))/velocitySquared;
          }
        }
        return time;
    }
    
    public void bump(AtomHard atom1, AtomHard atom2)
    {
        parentPhase.space.uEr1Mr2(r12, atom2.r, atom1.r);     //instance method      //r2-r1
        Space.uEa1Tv1Ma2Tv2(v12, atom2.rm, atom2.p, atom1.rm, atom1.p);  //v2-v1 = (p/m)2 - (p/m)1
        double factor = 2.0/(atom1.rm + atom2.rm)*Space.v1Dv2(r12,v12)/sig2;
        Space.uPEa1Tv1(atom1.p, factor, r12);
        Space.uMEa1Tv1(atom2.p, factor, r12);
    }
    
    public double getCollisionDiameter() {return collisionDiameter;}
    public void setCollisionDiameter(double c) {
        collisionDiameter = c;
        sig2 = c*c;
    }
}