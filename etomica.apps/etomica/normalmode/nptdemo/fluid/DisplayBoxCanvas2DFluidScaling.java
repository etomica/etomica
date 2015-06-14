package etomica.normalmode.nptdemo.fluid;

import java.awt.Graphics;

import etomica.action.activity.Controller;
import etomica.api.IAtom;
import etomica.api.IBox;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvas2D;
import etomica.space.ISpace;

public class DisplayBoxCanvas2DFluidScaling extends DisplayBoxCanvas2D {

    public DisplayBoxCanvas2DFluidScaling(DisplayBox _box, ISpace _space,
            Controller controller, IStuff stuff) {
        super(_box, _space, controller);
        p = _space.makeVector();
        this.stuff = stuff;
    }
    
    public void setPressure(double newPressure) {
        pressure = newPressure;
    }
    
    public double getPressure() {
        return pressure;
    }

    public void setDisplayDensity(double newDisplayDensity) {
        displayDensity = newDisplayDensity;
    }

    public double getDisplayDensity() {
        return displayDensity;
    }
    
    public void doPaint(Graphics g) {
        xScaling = stuff.stuff();
        super.doPaint(g);
    }
    
    protected void drawAtom(Graphics g, int[] origin, IAtom a) {

        IBox box = displayBox.getBox();
        double vOld = box.getBoundary().volume();
        int nAtoms = box.getLeafList().getAtomCount();
        double vNew = nAtoms/displayDensity;
        double rScale = Math.sqrt(vNew/vOld);
        
        p.Ea1Tv1((vNew-vOld)/rScale, xScaling[a.getLeafIndex()]);
        p.PE(a.getPosition());

        
        int sigmaP, xP, yP, baseXP, baseYP;

        g.setColor(displayBox.getColorScheme().getAtomColor(a));
        
        double toPixels = pixel.toPixels() * displayBox.getScale();

        baseXP = origin[0] + (int)(toPixels*p.getX(0));
        baseYP = origin[1] + (int)(toPixels*p.getX(1));
        /* Draw the core of the atom, specific to the dimension */
        double sigma = displayBox.getDiameterHash().getDiameter(a);
        // deafult diameter
        if (sigma == -1) sigma = 1;
        sigmaP = (int)(toPixels*sigma);
        sigmaP = (sigmaP == 0) ? 1 : sigmaP;
        xP = baseXP - (sigmaP>>1);
        yP = baseYP - (sigmaP>>1);
        g.fillOval(xP, yP, sigmaP, sigmaP);
    }

    protected double pressure, displayDensity;
    protected final IVectorMutable p;
    protected final IStuff stuff;
    protected IVector[] xScaling;
}
