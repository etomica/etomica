package simulate;

import java.beans.*;
import java.awt.*;

public class PotentialHardDisk extends Potential implements PotentialHard
{
    protected double collisionDiameter, sig2;

    public PotentialHardDisk(double d) {
        setCollisionDiameter(d);
    }

    public double collisionTime(AtomPair pair) {
        double r2 = pair.r2();
        double bij = pair.vDotr();
        if(r2 < sig2) {return (bij > 0) ? Double.MAX_VALUE : 0.0;}  //inside wall; no collision
        double time = Double.MAX_VALUE;

        if(bij < 0.0) {
          double velocitySquared = pair.v2();
          double discriminant = bij*bij - velocitySquared * ( r2 - sig2 );
          if(discriminant > 0) {
            time = (-bij - Math.sqrt(discriminant))/velocitySquared;
          }
        }
        return time;
    }
    
    public void bump(AtomPair pair) {
        double r2 = pair.r2();
        double factor = 2.0/(pair.atom1().rm() + pair.atom2().rm())*pair.vDotr()/r2;
        pair.cPair.push(factor);
    }
    
    public double getCollisionDiameter() {return collisionDiameter;}
    public void setCollisionDiameter(double c) {
        collisionDiameter = c;
        sig2 = c*c;
    }
    
    public double energy(AtomPair pair) {
        return (pair.r2() < sig2) ? Double.MAX_VALUE : 0.0;
    }
}