package simulate;
import java.util.Observable;
import java.util.Observer;

 /**
  * Container of information about the gravitational force acting on system.
  * Gravity is in effect only if an integrator and potential are used that look to it
  * in performing dynamics.
  * Changes in G can be propagated to interested objects by registering them as
  * observers of this class.
  * Default value of gravitational acceleration is zero.
  */

public class Gravity extends Observable {
    
    public static final double G_EARTH = 9.8*1e10/1e24;  //acceleration of gravity (on Earth), in A/ps^2
    private double g = 0.0;                  //acceleration of gravity for simulated system
    public double deltaG = 0.0;
    public double[] gVector = new Double[Space.D];

    public Gravity() {
        this(0.0);
    }
    
    public Gravity(double g) {
        this.g = g;
        Space.uEa1(gVector,0.0);
        gVector[Space.D-1] = g;
    }
    
    public void setG(double gNew) {
        deltaG = gNew - g;
        g = gNew;
        gVector[Space.D-1] = g;
        setChanged();
        notifyObservers();
    }
    public final double getG() {return g;}
}
    
