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

    public double collisionTime(AtomHard atom1, AtomHard atom2) {
   
        AtomHardDisk disk;
        AtomHardWall wall;
        if(atom2 instanceof AtomHardWall) {
           disk = (AtomHardDisk)atom1;
           wall = (AtomHardWall)atom2;
        }
        else {
           disk = (AtomHardDisk)atom2;
           wall = (AtomHardWall)atom1;
        }
        
        int i = (wall.isHorizontal()) ? 1 : 0;  //indicates if collision affects x or y coordinate
        
        //wall or disk has non-zero force
        if(!wall.isForceFree() || !disk.isForceFree()) {  
            return timeWithAcceleration(i, disk, wall);
        }
        //wall or disk is stationary and gravity is acting       
        else if(!parentPhase.noGravity && i==0 && (wall.isStationary() || disk.isStationary())) {  
            return timeWithAcceleration(i, disk, wall);
        }
        else {return timeNoAcceleration(i, disk, wall);}
    }
    
    private double timeNoAcceleration(int i, AtomHardDisk disk, AtomHardWall wall) {
      double dr = parentPhase.space.r1iMr2i(i,wall.r,disk.r);
      double dv = wall.p[i]*wall.rm-disk.p[i]*disk.rm;
       
      if(dr*dv > 0.0) {return Double.MAX_VALUE;}    //Separating, no collision

      double adr = Math.abs(dr);
      if(adr < collisionRadius) {return 0.0;}            // Inside core and approaching; collision now
      else {return (adr-collisionRadius)/Math.abs(dv);}  //Outside core and approaching; collision at core
    }
    
    private double timeWithAcceleration(int i, AtomHardDisk disk, AtomHardWall wall) {
      double time = Double.MAX_VALUE;
      double dr, dv;
      dr = parentPhase.space.r1iMr2i(i,wall.r,disk.r);
      dv = wall.p[i]*wall.rm - disk.p[i]*disk.rm;
      
      if(Math.abs(dr) < collisionRadius) {   //inside wall; collision now if approaching, never otherwise
        return (dr*dv > 0) ? Double.MAX_VALUE : 0.0;}
        
      dr += (dr > 0.0) ? -collisionRadius : +collisionRadius;
      double a = wall.f[i]*wall.rm - disk.f[i]*disk.rm;
      if(i==1) {  //consider acceleration of gravity, which acts on non-stationary atom and/or wall
        if(!wall.isStationary()) {a += parentPhase.getG();}  
        if(!disk.isStationary()) {a -= parentPhase.getG();}
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
    
    public void bump(AtomHard atom1, AtomHard atom2)  //this needs updating to check for isStationary
    {
        AtomHardDisk disk;
        AtomHardWall wall;
        if(atom2 instanceof AtomHardWall) {
           disk = (AtomHardDisk)atom1;
           wall = (AtomHardWall)atom2;
        }
        else {
           disk = (AtomHardDisk)atom2;
           wall = (AtomHardWall)atom1;
        }
        
        int i = (wall.isHorizontal()) ? 1 : 0;  //indicates if collision affects x or y coordinate

        if(wall.isStationary()) {
            disk.p[i] *= -1;
            wall.accumulateP(2*Math.abs(disk.p[i]));
        }
        else {
          double dv = wall.p[i]*wall.rm - disk.p[i]*disk.rm;
          double dp = -2.0/(wall.rm + disk.rm)*dv;
          wall.p[i] += dp;
          disk.p[i] -= dp; 
        }
    }
    
    public double getCollisionDiameter() {return collisionDiameter;}
    public void setCollisionDiameter(double c) {
        collisionDiameter = c;
        collisionRadius = 0.5*c;
    }
}