package simulate;
import java.util.Observable;
import java.util.Observer;
import java.util.*;

 /**
  * Container of information about the gravitational force acting on system.
  * Gravity is in effect only if an integrator and potential are used that look to it
  * in performing dynamics.
  * Changes in G can be propagated to interested objects by registering them as
  * observers of this class.
  * Default value of gravitational acceleration is zero.
  */

public class Gravity extends java.util.Observable      {
    
    public static final double G_EARTH = 9.8*1e10/1e24;  //acceleration of gravity (on Earth), in A/ps^2
    private double g = 0.0;                  //acceleration of gravity for simulated system
    public double deltaG = 0.0;
    public Space.Vector gVector;

    public Gravity(double g) {
        this.g = g;
        gVector = new Space2D.Vector();  //assume 2D space for now
        gVector.E(0.0);
        gVector.setComponent(Simulation.D-1,g);
    }
    
    public void setG(double gNew) {
        deltaG = gNew - g;
        g = gNew;
        gVector.setComponent(Simulation.D-1,g);
        setChanged();
        notifyObservers();
    }
    public final double getG() {return g;}
}
    
