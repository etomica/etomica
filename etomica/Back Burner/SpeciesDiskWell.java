package simulate;
import java.util.*;
import java.awt.*;
import java.beans.*;

public class SpeciesDiskWell extends SpeciesDisk {
  double lambda = 1.5;
  Color wellColor = Color.gray;

  public SpeciesDiskWell() {
    super();
  }

  public void draw(Graphics g, int[] origin, double scale) {
    double toPixels = scale*Phase.TO_PIXELS;
    int sigmaP = (int)(toPixels*sigma*lambda);
    g.setColor(wellColor);
    double radius = halfSigma*lambda;
    for(int i=nElements; --i>=0; ) {
        int xP = origin[0] + (int)(toPixels*(molecule[i].r[0]-radius));
        int yP = origin[1] + (int)(toPixels*(molecule[i].r[1]-radius));
        g.drawOval(xP,yP,sigmaP,sigmaP);
    }
    sigmaP = (int)(toPixels*sigma);
    g.setColor(color);
    for(int i=nElements; --i>=0; ) {
        int xP = origin[0] + (int)(toPixels*(molecule[i].r[0]-halfSigma));
        int yP = origin[1] + (int)(toPixels*(molecule[i].r[1]-halfSigma));
        g.fillOval(xP,yP,sigmaP,sigmaP);
    }
  }

  public void setLambda(double lam) {lambda = lam;}
  public double getLambda() {return lambda;}

  public void setWellColor(Color c) {wellColor = c;}
  public Color getWellColor() {return wellColor;}

}