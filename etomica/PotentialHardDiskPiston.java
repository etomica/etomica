package simulate;

import simulate.*;
import java.beans.*;
import java.awt.*;

public class PotentialHardDiskPiston extends simulate.Potential
{
    private double collisionDiameter;
    private double collisionRadius;

    public PotentialHardDiskPiston(double d) {
        setCollisionDiameter(d);
    }

    public double collisionTime(Atom atom1, Atom atom2) {
   
        Atom disk;
        AtomWall wall;
        if(atom2 instanceof AtomWall) {
           disk = atom1;
           wall = (AtomWall)atom2;
        }
        else {
           disk = atom2;
           wall = (AtomWall)atom1;
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
*/            return time;}  //inside wall; no collision
        dr += (dr > 0.0) ? -collisionRadius : +collisionRadius;
        dv = wall.p[i]*wall.rm - disk.p[i]*disk.rm;
        double a = wall.f[i]*wall.rm;
        double discrim = dv*dv - 2*a*dr;
        if(discrim > 0) {
            boolean adr = (a*dr > 0);
            boolean adv = (a*dv > 0);
            if(adr && adv) {time = Double.MAX_VALUE;}
            else if(adr) {time = (-dv - Math.sqrt(discrim))/a;}
            else if(-a*dr/(dv*dv) < 1.e-7) {if(dr*dv<0) time = -dr/dv*(1+0.5*dr*a/(dv*dv));}
            else {time = (-dv + Math.sqrt(discrim))/a;}
        }
//        System.out.println(dr+"  "+dv+"  "+a+"  "+time);
        return time;
    }
    
    public void bump(Atom atom1, Atom atom2)
    {
        Atom disk;
        AtomWall wall;
        if(atom2 instanceof AtomWall) {
           disk = atom1;
           wall = (AtomWall)atom2;
        }
        else {
           disk = atom2;
           wall = (AtomWall)atom1;
        }

        int i;
        if(wall.isVertical()) {i = 0;}
        else if(wall.isHorizontal()) {i = 1;}
        else {i = 0;}

        double dv = wall.p[i]*wall.rm - disk.p[i]*disk.rm;
        double dp = -2.0/(atom1.rm + atom2.rm)*dv;
        wall.p[i] += dp;
        disk.p[i] -= dp;        
    }
    
    public double getCollisionDiameter() {return collisionDiameter;}
    public void setCollisionDiameter(double c) {
        collisionDiameter = c;
        collisionRadius = 0.5*c;
    }
}