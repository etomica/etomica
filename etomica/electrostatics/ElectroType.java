package etomica.electrostatics;

import etomica.units.*;
import etomica.Space;
import java.awt.Graphics;

/**
 * Parent class for all electrostatic atom types (monopole, point multipoles)
 */
public abstract class ElectroType implements java.io.Serializable {
    
    //NOTE:  Subclasses currently are implementing draw incorrectly.
    // They are written using scale as the double parameter, whereas in the atom's draw
    // method it is passing the length toPixels value, which includes scale.
    // electrotypes use scale and the TO_PIXELS from the corresponding unit class 
    public abstract void draw(Graphics g, int origin[], double t, Space.Vector r);
    public static class Null extends ElectroType {
        public void draw(Graphics g, int origin[], double t, Space.Vector r) {}
    }
    
}