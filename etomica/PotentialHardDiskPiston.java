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
        
        double time = Double.MAX_VALUE;
        double dr, dv;
        int i;
        
        AtomType.Wall wallType = (AtomType.Wall)wall.type;
        if(wallType.isVertical()) {i = 0;}
        else {i = 1;}
        
        dr = wall.position(i) - disk.position(i);
        if(Math.abs(dr) < collisionRadius) { 
 /*           double eps = 1.e-5;
            wall.r[i] = 0.5*((2+eps)*wall.r[i] - eps*disk.r[i]);
            disk.r[i] = 0.5*((2+eps)*disk.r[i] - eps*wall.r[i]);
            dr *= 1+eps;
        }
*/            return 0.0;}  //inside wall; no collision
        dr += (dr > 0.0) ? -collisionRadius : +collisionRadius;
        dv = wall.momentum(i)*wall.rm() - disk.momentum(i)*disk.rm();
//        double a = wall.f[i]*wall.rm;  //acceleration of wall
        double a = 0.0;   //temporary
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
    
    public void bump(AtomPair pair)
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

        int i;
        AtomType.Wall wallType = (AtomType.Wall)wall.type;
        if(wallType.isVertical()) {i = 0;}
        else if(wallType.isHorizontal()) {i = 1;}
        else {i = 0;}

        double dv = wall.momentum(i)*wall.rm() - disk.momentum(i)*disk.rm();
        double dp = -2.0/(pair.atom1().rm() + pair.atom2().rm())*dv;
        if(!wall.isStationary()) wall.momentum().PE(i,dp);
        if(!disk.isStationary()) disk.momentum().PE(i,-dp);
        double dr = wall.position(i) - disk.position(i);
//        double dr = parentPhase.space.r1iMr2i(i,wall.r,disk.r);
//        disk.r[i] = wall.r[i] - (1.+1.e-5)*collisionRadius*Math.abs(dr)/dr;
        disk.position().PE(i,-1.e-5*dr);
    }
    
    public double getCollisionDiameter() {return collisionDiameter;}
    public void setCollisionDiameter(double c) {
        collisionDiameter = c;
        collisionRadius = 0.5*c;
    }
}