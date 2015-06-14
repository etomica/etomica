package etomica.normalmode.nptdemo;

import java.awt.Graphics;

import etomica.action.activity.Controller;
import etomica.api.IAtom;
import etomica.api.IBox;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvas2D;
import etomica.normalmode.CoordinateDefinition;
import etomica.space.ISpace;

public class DisplayBoxCanvas2DNpTScaling extends DisplayBoxCanvas2D {

    public DisplayBoxCanvas2DNpTScaling(DisplayBox _box, ISpace _space,
            Controller controller, CoordinateDefinition coordinateDefinition) {
        super(_box, _space, controller);
        p = _space.makeVector();
        this.coordinateDefinition = coordinateDefinition;
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
    
    protected void drawAtom(Graphics g, int[] origin, IAtom a) {
        p.E(a.getPosition());
        
        IVector l = coordinateDefinition.getLatticePosition(a);
        p.ME(l);

        IBox box = displayBox.getBox();
        double vOld = box.getBoundary().volume();
        int nAtoms = box.getLeafList().getAtomCount();
        double vNew = nAtoms/displayDensity;
        double rScale = Math.sqrt(vNew/vOld);
        double latticeScale = Math.exp((pressure*(vNew-vOld))/((nAtoms-1)*1*2))/rScale;
        
        p.TE(latticeScale);
        p.PE(l);
        
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
    protected final CoordinateDefinition coordinateDefinition;
}
