package simulate;

import java.beans.*;
import java.awt.*;

// Written only for stationary wall

public class PotentialHardDiskWall extends Potential implements PotentialHard
{
    protected double collisionDiameter, collisionRadius;

    public PotentialHardDiskWall(double d) {
        setCollisionDiameter(d);
    }

    public double collisionTime(AtomPair pair) {
   
        Atom disk;
        Atom wall;
        if(pair.atom2().type instanceof AtomType.Wall) {
           disk = pair.atom1();
           wall = pair.atom2();
        }
        else {
           disk = pair.atom2;
           wall = pair.atom1;
        }
        AtomType.Wall wallType = (AtomType.Wall)wall.type;
                
        int i = (((AtomType.Wall)wall.type).isHorizontal()) ? 1 : 0;  //indicates if collision affects x or y coordinate
        
        //wall or disk has non-zero force
 //       if(!wall.isForceFree() || !disk.isForceFree()) {  
 //           return timeWithAcceleration(i, disk, wall);
 //       }
        //wall or disk is stationary and gravity is acting       
 //       else if(!parentPhase.noGravity && i==0 && (wall.isStationary() || disk.isStationary())) {  
 //           return timeWithAcceleration(i, disk, wall);
 //       }
 //       else {
            return timeNoAcceleration(i, disk, wall);//}
    }
    
    private double timeNoAcceleration(int i, Atom disk, Atom wall) {
//      double dr = parentPhase.space.r1iMr2i(i,wall.r,disk.r);
      double dr = wall.position(i) - disk.position(i);   //no PBC
      double dv = wall.momentum(i)*wall.rm()-disk.momentum(i)*disk.rm();
       
      if(dr*dv > 0.0) {return Double.MAX_VALUE;}    //Separating, no collision

      double adr = Math.abs(dr);
      if(adr < collisionRadius) {return 0.0;}            // Inside core and approaching; collision now
      else {return (adr-collisionRadius)/Math.abs(dv);}  //Outside core and approaching; collision at core
    }
    
    private double timeWithAcceleration(int i, Atom disk, Atom wall) {
      double time = Double.MAX_VALUE;
//      double dr, dv;
//      dr = parentPhase.space.r1iMr2i(i,wall.r,disk.r);
//      dv = wall.p[i]*wall.rm - disk.p[i]*disk.rm;
      double dr = wall.position(i) - disk.position(i);   //no PBC
      double dv = wall.momentum(i)*wall.rm()-disk.momentum(i)*disk.rm();
      
      if(Math.abs(dr) < collisionRadius) {   //inside wall; collision now if approaching, never otherwise
        return (dr*dv > 0) ? Double.MAX_VALUE : 0.0;}
        
      dr += (dr > 0.0) ? -collisionRadius : +collisionRadius;
//      double a = wall.f[i]*wall.rm - disk.f[i]*disk.rm;
      double a = 0.0;  //this needs to be changed to the line above
      if(i==1) {  //consider acceleration of gravity, which acts on non-stationary atom and/or wall
        if(!wall.isStationary()) {a += wall.parentPhase().getG();}  
        if(!disk.isStationary()) {a -= disk.parentPhase().getG();}
      }
      double discrim = dv*dv - 2*a*dr;
      if(discrim > 0) {
        boolean adr = (a*dr > 0);
        boolean adv = (a*dv > 0);
        int aSign = (a > 0) ? +1 : -1;
        if(adr && adv) {time = Double.MAX_VALUE;}
        else if(adr) {time = (-dv - aSign*Math.sqrt(discrim))/a;}
        else if(-a*dr/(dv*dv) < 1.e-7) {if(dr*dv<0) time = -dr/dv*(1+0.5*dr*a/(dv*dv));} //series expansion for small acceleration
        else {time = (-dv + aSign*Math.sqrt(discrim))/a;}
      }
      return time;
    }
    
    public void bump(AtomPair pair)  //this needs updating to check for isStationary
    {
        Atom disk;
        Atom wall;
        if(pair.atom2().type instanceof AtomType.Wall) {
           disk = pair.atom1();
           wall = pair.atom2();
        }
        else {
           disk = pair.atom2;
           wall = pair.atom1;
        }
        AtomType.Wall wallType = (AtomType.Wall)wall.type;
                
        int i = (((AtomType.Wall)wall.type).isHorizontal()) ? 1 : 0;  //indicates if collision affects x or y coordinate

        if(wall.isStationary()) {
            disk.momentum().TE(i,-1.0);
//            wall.accumulateP(2*Math.abs(disk.momentum(i)));
        }
        else {
          double dv = wall.momentum(i)*wall.rm()-disk.momentum(i)*disk.rm();
          double dp = -2.0/(wall.rm() + disk.rm())*dv;
          wall.momentum().PE(i,+dp);
          disk.momentum().PE(i,-dp); 
        }
    }
    
    public double getCollisionDiameter() {return collisionDiameter;}
    public void setCollisionDiameter(double c) {
        collisionDiameter = c;
        collisionRadius = 0.5*c;
    }
}