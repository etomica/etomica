package simulate;

import java.beans.*;
import java.awt.*;

// Written only for stationary wall

public class PotentialHardDiskWall extends Potential implements PotentialHard
{
    double collisionDiameter, collisionRadius, sig2;

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
        
        double time = Double.MAX_VALUE;
        int i;
        
        if(wall.isVertical()) {i = 0;}
        else if(wall.isHorizontal()) {i = 1;}
        else {i = 0;}
        
        if(parentPhase.noGravity || i==0 || !(wall.isStationary() || disk.isStationary())) {
            double dr, t, dtdr;
            dr = parentPhase.space.r1iMr2i(i,wall.r,disk.r);
            dtdr = 1.0/(disk.p[i]*disk.rm);
            t = dr*dtdr;

            if(t > 0.0) {time = Math.max(0.0,t-collisionDiameter/2.0*Math.abs(dtdr));}
            return time;
        }
        else {
            double dr, dv;
            dr = parentPhase.space.r1iMr2i(i,wall.r,disk.r);
            dv = wall.p[i]*wall.rm - disk.p[i]*disk.rm;
            if(Math.abs(dr) < collisionRadius) {   //this may still need some work
                return (dr*dv > 0) ? Double.MAX_VALUE : 0.0;}  //inside wall; no collision
            dr += (dr > 0.0) ? -collisionRadius : +collisionRadius;
            double a = wall.isStationary() ? -(parentPhase.getG()) : parentPhase.getG();
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
        
        if(wall.isVertical()) {
            disk.p[0] *= -1;
            wall.accumulateP(2*Math.abs(disk.p[0]));
        }
        if(wall.isHorizontal()) {
            disk.p[1] *= -1;
            wall.accumulateP(2*Math.abs(disk.p[1]));
        }
        
    }
    
    public double getCollisionDiameter() {return collisionDiameter;}
    public void setCollisionDiameter(double c) {
        collisionDiameter = c;
        collisionRadius = 0.5*c;
        sig2 = c*c;
    }
}