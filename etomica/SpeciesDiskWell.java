package simulate;
import java.util.*;
import java.awt.*;
import java.beans.*;

public class SpeciesDiskWell extends SpeciesDisks {
  double lambda;
  Color wellColor;

  public SpeciesDiskWell() {
    super();
  }
  
  public void setDefaults() {
    super.setDefaults();
    setWellColor(Color.gray);
    setLambda(1.5);
  }

  void initializeMolecules() {
    initializeMolecules(diameter, mass, color, lambda, wellColor);
  }
  void initializeMolecules(double d, double m, Color c, double l, Color w) {
    setDiameter(d);    //call set methods to pass diameter and mass to atoms
    setMass(m);
    setColor(c);
    setLambda(l);
    setWellColor(w);
  }

  public void draw(Graphics g, int[] origin, double scale) {
    double toPixels = scale*Phase.TO_PIXELS;
    double halfWell = radius*lambda;
    Atom nextSpeciesAtom = lastAtom.getNextAtom();
    int wellP = (int)(toPixels*diameter*lambda);
    int diameterP = (int)(toPixels*diameter);
    
    g.setColor(wellColor);
    for(Atom a=firstAtom; a!=nextSpeciesAtom; a=a.getNextAtom()) {
        int xP = origin[0] + (int)(toPixels*(a.r[0]-halfWell));
        int yP = origin[1] + (int)(toPixels*(a.r[1]-halfWell));
        g.drawOval(xP,yP,wellP,wellP);
    }

    if(parentPhase.drawOverflowImages) {
        for(Atom a=firstAtom; a!=nextSpeciesAtom; a=a.getNextAtom()) {
            double[][] shifts = parentPhase.space.getOverflowShifts(a.r,halfWell);
            for(int i=0; i<shifts.length; i++) {
               int xP = origin[0] + (int)(toPixels*(shifts[i][0]+a.r[0]-halfWell));
               int yP = origin[1] + (int)(toPixels*(shifts[i][1]+a.r[1]-halfWell));
               g.setColor(wellColor);
               g.drawOval(xP,yP,wellP,wellP);
               xP = origin[0] + (int)(toPixels*(shifts[i][0]+a.r[0]-radius));
               yP = origin[1] + (int)(toPixels*(shifts[i][1]+a.r[1]-radius));
               g.setColor(color);
               g.fillOval(xP,yP,diameterP,diameterP);
            }
        }
    }
    g.setColor(color);
    for(Atom a=firstAtom; a!=nextSpeciesAtom; a=a.getNextAtom()) {
        int xP = origin[0] + (int)(toPixels*(a.r[0]-radius));
        int yP = origin[1] + (int)(toPixels*(a.r[1]-radius));
        g.fillOval(xP,yP,diameterP,diameterP);
    }         
  }

  public double getLambda() {return lambda;}
  public void setLambda(double lam) {lambda = lam;}

  public Color getWellColor() {return wellColor;}
  public void setWellColor(Color c) {wellColor = c;}
}