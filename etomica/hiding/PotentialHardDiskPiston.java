package simulate;

import java.beans.*;
import java.awt.*;

public class PotentialHardDiskPiston extends Potential implements PotentialHard
{
    private double collisionDiameter;
    private double collisionRadius;

    public PotentialHardDiskPiston(double d) {
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
        double dr, dv;
        int i;
        
        if(wall.isVertical()) {i = 0;}
        else if(wall.isHorizontal()) {i = 1;}
        else {i = 0;}
        
        dr = wall.r[i] - disk.r[i];
        if(Math.abs(dr) < collisionRadius) { 
 /*           double eps = 1.e-5;
            wall.r[i] = 0.5*((2+eps)*wall.r[i] - eps*disk.r[i]);
            disk.r[i] = 0.5*((2+eps)*disk.r[i] - eps*wall.r[i]);
            dr *= 1+eps;
        }
*/            return 0.0;}  //inside wall; no collision
        dr += (dr > 0.0) ? -collisionRadius : +collisionRadius;
        dv = wall.p[i]*wall.rm - disk.p[i]*disk.rm;
        double a = wall.f[i]*wall.rm;
        double discrim = dv*dv - 2*a*dr;
        if(discrim > 0) {
            boolean adr = (a*dr > 0);
            boolean adv = (a*dv > 0);
            int aSign = (a > 0) ? +1 : -1;
            if(adr && adv) {time = Double.MAX_VALUE;}
            else if(adr) {time = (-dv - aSign*Math.sqrt(discrim))/a;}
            else if(-a*dr/(dv*dv) < 1.e-7) {if(dr*dv<0) time = -dr/dv*(1+0.5*dr*a/(dv*dv));}
            else {time = (-dv + aSign*Math.sqrt(discrim))/a;}
        }
//        System.out.println(dr+"  "+dv+"  "+a+"  "+time);
        return time;
    }
    
    public void bump(AtomHard atom1, AtomHard atom2)
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

        int i;
        if(wall.isVertical()) {i = 0;}
        else if(wall.isHorizontal()) {i = 1;}
        else {i = 0;}

        double dv = wall.p[i]*wall.rm - disk.p[i]*disk.rm;
        double dp = -2.0/(atom1.rm + atom2.rm)*dv;
        wall.p[i] += dp;
        disk.p[i] -= dp; 
        double dr = parentPhase.space.r1iMr2i(i,wall.r,disk.r);
        disk.r[i] = wall.r[i] - (1.+1.e-5)*collisionRadius*Math.abs(dr)/dr;
    }
    
    public double getCollisionDiameter() {return collisionDiameter;}
    public void setCollisionDiameter(double c) {
        collisionDiameter = c;
        collisionRadius = 0.5*c;
    }
}