package etomica.electrostatics;

import etomica.units.*;
import etomica.Space;
import java.awt.Graphics;

    /**
     * ElectroType for a simple Coulombic charge (monopole) 
     * Default units are electron charge
     */
public final class Monopole extends ElectroType {
    private double z;  //coulomb charge
            
    public Monopole() {
        setZ(0.0);
    }
    public Monopole(double z) {
        setZ(z);
    }
    public void setZ(double t) {z = t;}
    public final double getZ() {return z;}
    public Dimension getZDimension() {return Dimension.CHARGE;}
    
    public void draw(Graphics g, int origin[], double scale, Space.Vector r) {}
}
