package simulate;
import java.util.*;
import java.awt.*;
import java.beans.*;

public class SpeciesDiskWell extends SpeciesDisks {
  double lambda;
  Color wellColor;

  public SpeciesDiskWell(PhaseSpace ps) {
    this(ps, 20,1);
  }
  public SpeciesDiskWell(PhaseSpace ps, int nM, int nA) {
    super(/*ps, */nM, nA, new AtomType.Well(1.0,Color.black,0.1,1.5));
  }
  
/*  public void draw(Graphics g, int[] origin, double scale) {
    double toPixels = scale*DisplayConfiguration.SIM2PIXELS;
    double halfWell = radius*lambda;
    Atom nextSpeciesAtom = lastAtom().nextAtom();
    int wellP = (int)(toPixels*diameter*lambda);
    int diameterP = (int)(toPixels*diameter);
    
    g.setColor(wellColor);
    for(AtomC a=(AtomC)firstAtom(); a!=nextSpeciesAtom; a=(AtomC)a.getNextAtom()) {
        int xP = origin[0] + (int)(toPixels*(a.r[0]-halfWell));
        int yP = origin[1] + (int)(toPixels*(a.r[1]-halfWell));
        g.drawOval(xP,yP,wellP,wellP);
    }

    if(DisplayConfiguration.DRAW_OVERFLOW) {
        for(AtomC a=(AtomC)firstAtom(); a!=nextSpeciesAtom; a=(AtomC)a.getNextAtom()) {
            double[][] shifts = parentSimulation.space.getOverflowShifts(a.r,halfWell);
            for(int i=0; i<shifts.length; i++) {
               int xP = origin[0] + (int)(toPixels*(shifts[i][0]+a.r[0]-halfWell));
               int yP = origin[1] + (int)(toPixels*(shifts[i][1]+a.r[1]-halfWell));
               g.setColor(wellColor);
               g.drawOval(xP,yP,wellP,wellP);
               xP = origin[0] + (int)(toPixels*(shifts[i][0]+a.r[0]-radius));
               yP = origin[1] + (int)(toPixels*(shifts[i][1]+a.r[1]-radius));
               colorScheme.setAtomColor(a);
               g.setColor(a.getColor());
               g.fillOval(xP,yP,diameterP,diameterP);
            }
        }
    }
    for(AtomC a=(AtomC)firstAtom(); a!=nextSpeciesAtom; a=(AtomC)a.getNextAtom()) {
        int xP = origin[0] + (int)(toPixels*(a.r[0]-radius));
        int yP = origin[1] + (int)(toPixels*(a.r[1]-radius));
        colorScheme.setAtomColor(a);
        g.setColor(a.getColor());
        g.fillOval(xP,yP,diameterP,diameterP);
    }         
  }*/

  public double getLambda() {return ((AtomType.Well)protoType).lambda();}
  public void setLambda(double lam) {((AtomType.Well)protoType).setLambda(lam);}
}