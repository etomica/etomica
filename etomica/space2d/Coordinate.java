package simulate.space2D;
import simulate.Space;

public class Coordinate extends Space.Coordinate {
    public final Vector r = new Vector();  //Cartesian coordinates
    public final Vector p = new Vector();  //Momentum vector
    public Coordinate(Space.Occupant o) {super(o);}
    public Space.Vector position() {return r;}
    public Space.Vector momentum() {return p;}
    public double position(int i) {return r.component(i);}
    public double momentum(int i) {return p.component(i);}
    public final double kineticEnergy(double mass) {return 0.5*p.squared()/mass;}
} 
