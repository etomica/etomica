package simulate.electrostatics;

import simulate.units.*;
import simulate.*;
import java.awt.Graphics;
import java.awt.Color;

//draw method is specific to 2D

    /**
     * ElectroType for a point dipole 
     * Default units are Debyes
     */
public final class Dipole extends ElectroType {
    private double mu;  //magnitude
    private final Space.Vector e; //direction
            
    public Dipole() {
        e = Simulation.space().makeVector();
        initializeE();
        setMu(0.0);
    }
    public Dipole(double mu) {
        e = Simulation.space().makeVector();
        initializeE();
        setMu(mu);
    }
    
    private void initializeE() {
        e.E(0.0);   //zero vector
        e.setComponent(0,1.0); //set first component to unity
    }
    
    public Space.Vector e() {return e;}
    
    public void setMu(double t) {mu = t;}
    public final double getMu() {return mu;}
    public Dimension getMuDimension() {return Dimension.DIPOLE;}
    
    
    public void draw(Graphics g, int origin[], double scale, Space.Vector r) {
        double length = scale * BaseUnit.Dipole.Sim.UNIT.toPixels(mu);
        int circleRadius = (int)(0.2*length);
        int x0 = origin[0] + (int)BaseUnit.Length.Sim.UNIT.toPixels(scale*r.component(0));
        int y0 = origin[1] + (int)BaseUnit.Length.Sim.UNIT.toPixels(scale*r.component(1));
        int xa = (int)(x0 - 0.5*length*e.component(0));
        int ya = (int)(y0 - 0.5*length*e.component(1));
        int xb = (int)(x0 + 0.5*length*e.component(0));
        int yb = (int)(y0 + 0.5*length*e.component(1));
        g.setColor(Color.black);
        g.drawLine(xa, ya, xb, yb);
        g.setColor(Color.green);
        xb -= circleRadius*(1+e.component(0));
        yb -= circleRadius*(1+e.component(1));
        g.fillOval(xb, yb, 2*circleRadius, 2*circleRadius);
    }
}
