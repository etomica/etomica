package simulate.space2D;
import simulate.Space;
import simulate.Phase;
import java.util.Random;
import java.awt.Graphics;

   /**
    * Class for implementing no periodic boundary conditions
    */
public class BoundaryNone extends Boundary {
    private final Vector temp = new Vector();
    private final double[][] shift0 = new double[0][Space2D.D];
    public final Vector dimensions = new Vector();
    public Space.Vector dimensions() {return dimensions;}
    public static final Random random = new Random();
    public BoundaryNone() {super();}
    public BoundaryNone(Phase p) {super(p);}
    public final void nearestImage(Space.Vector dr) {}
    public final void centralImage(Space.Vector r) {}
    public final void nearestImage(Vector dr) {}
    public final void centralImage(Vector r) {}
    public final void centralImage(Coordinate c) {}
    public double volume() {return Double.MAX_VALUE;}
    public void inflate(double s) {}
    public double[][] imageOrigins(int nShells) {return new double[0][Space2D.D];}
    public double[][] getOverflowShifts(Space.Vector rr, double distance) {return shift0;}
    public Space.Vector randomPosition() {  //arbitrary choice for this method in this boundary
        temp.x = random.nextDouble(); 
        temp.y = random.nextDouble(); 
        return temp;
    }
    public void draw(Graphics g, int[] origin, double scale) {}
}
    
