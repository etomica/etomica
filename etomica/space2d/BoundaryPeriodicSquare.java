package simulate.space2D;
import simulate.Space;
import simulate.Phase;
import simulate.Potential;
import simulate.AtomPair;
import simulate.Atom;
import simulate.AtomType;
import java.awt.Graphics;
import java.awt.Color;
import simulate.units.BaseUnit;
import java.util.Random;
import simulate.Default;

   /**
    * Class for implementing rectangular periodic boundary conditions
    */
public class BoundaryPeriodicSquare extends Boundary implements Potential.Hard {
    private final Vector temp = new Vector();
    public static final Random random = new Random();
    private final double[][] shift0 = new double[0][Space2D.D];
    private final double[][] shift1 = new double[1][Space2D.D]; //used by getOverflowShifts
    private final double[][] shift3 = new double[3][Space2D.D];
    public BoundaryPeriodicSquare() {this(Default.BOX_SIZE,Default.BOX_SIZE);}
    public BoundaryPeriodicSquare(Phase p) {super(p);}
    public BoundaryPeriodicSquare(Phase p, double lx, double ly) {super(p); dimensions.x = lx; dimensions.y = ly;}
    public BoundaryPeriodicSquare(double lx, double ly) {dimensions.x = lx; dimensions.y = ly;}
    public final Vector dimensions = new Vector();
    public final Space.Vector dimensions() {return dimensions;}
    public Space.Vector randomPosition() {
        temp.x = dimensions.x*random.nextDouble(); 
        temp.y = dimensions.y*random.nextDouble(); 
        return temp;}
    public void nearestImage(Space.Vector dr) {nearestImage((Vector)dr);}
    public void nearestImage(Vector dr) {
        dr.x -= dimensions.x * ((dr.x > 0.0) ? Math.floor(dr.x/dimensions.x+0.5) : Math.ceil(dr.x/dimensions.x-0.5));
        dr.y -= dimensions.y * ((dr.y > 0.0) ? Math.floor(dr.y/dimensions.y+0.5) : Math.ceil(dr.y/dimensions.y-0.5));
    }
    public void centralImage(Coordinate c) {centralImage(c.r);}
    public void centralImage(Space.Vector r) {centralImage((Vector)r);}
    public void centralImage(Vector r) {
        r.x -= dimensions.x * ((r.x >= 0.0) ? Math.floor(r.x/dimensions.x) : Math.ceil(r.x/dimensions.x-1.0));
        r.y -= dimensions.y * ((r.y >= 0.0) ? Math.floor(r.y/dimensions.y) : Math.ceil(r.y/dimensions.y-1.0));
    }
    public void inflate(double scale) {dimensions.TE(scale);}
    public double volume() {return dimensions.x * dimensions.y;}
        
    //"Collision" whenever atom travels half the edge length of the simulation volume
    //needs some work to handle non-disk atoms better
    public double collisionTime(AtomPair pair) {
        Atom a = pair.atom1();
        if(!(a.type instanceof AtomType.Disk)) {
            return Double.MAX_VALUE;}
        Vector p = (Vector)a.coordinate().momentum();
        double diameter = ((AtomType.Disk)a.type).diameter();
        //assumes range of potential is .le. diameter, simulation box is square (or x is smaller dimension)
        return 0.5*(dimensions.y-1.0001*diameter)/(a.rm()*Math.sqrt(p.squared()));  
    }
    //No action needed at collision, just want to update neighbor list
    public void bump(AtomPair pair) {}
    public double lastCollisionVirial() {return 0.0;}
    public Space.Tensor lastCollisionVirialTensor() {return Tensor.ZERO;}
    public double energy(AtomPair pair) {return 0.0;}
    public double energyLRC(int n1, int n2, double V) {return 0.0;}
    public boolean overlap(AtomPair pair) {return false;}

    public void draw(Graphics g, int[] origin, double scale) {
        g.setColor(Color.gray);
        double toPixels = scale*BaseUnit.Length.Sim.TO_PIXELS;
        g.drawRect(origin[0],origin[1],(int)(toPixels*dimensions.component(0))-1,(int)(toPixels*dimensions.component(1))-1);
        }
    /** Computes origins for periodic images
    */
    public double[][] imageOrigins(int nShells) {
        int nImages = (2*nShells+1)*(2*nShells+1)-1;
        double[][] origins = new double[nImages][Space2D.D];
        int k = 0;
        for(int i=-nShells; i<=nShells; i++) {
            for(int j=-nShells; j<=nShells; j++) {
                if(i==0 && j==0) {continue;}
                origins[k][0] = i*dimensions.x;
                origins[k][1] = j*dimensions.y;
                k++;
            }
        }
        return origins;
    }

    /** Returns coordinate shifts needed to draw all images that overflow into central image
        * 0, 1, or 3 shifts may be returned
        */
    public double[][] getOverflowShifts(Space.Vector rr, double distance) {
        Vector r = (Vector)rr;
        int shiftX = 0;
        int shiftY = 0;
        if(r.x-distance < 0.0) {shiftX = +1;}
        else if(r.x+distance > dimensions.x) {shiftX = -1;}
            
        if(r.y-distance < 0.0) {shiftY = +1;}
        else if(r.y+distance > dimensions.y) {shiftY = -1;}
            
        if(shiftX == 0) {
            if(shiftY == 0) {
                return shift0;
            }
            else {
                shift1[0][0] = 0.0;
                shift1[0][1] = shiftY*dimensions.y;
                return shift1;
            }
        }
        else { //shiftX != 0
            if(shiftY == 0) {
                shift1[0][0] = shiftX*dimensions.x;
                shift1[0][1] = 0.0;
                return shift1;
            }
            else {
                shift3[0][0] = shiftX*dimensions.x;
                shift3[0][1] = 0.0;
                shift3[1][0] = 0.0;
                shift3[1][1] = shiftY*dimensions.y;
                shift3[2][0] = shift3[0][0];
                shift3[2][1] = shift3[1][1];
                return shift3;
            }
        }
    } //end of getOverflowShifts
}  //end of BoundaryPeriodicSquare
