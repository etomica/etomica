package simulate.space2D;
import simulate.Space;
import simulate.Atom;
import simulate.Phase;
import simulate.AtomPair;

public final class BoundaryHard extends BoundaryPeriodicSquare {
    private double collisionRadius = 0.0;
    public BoundaryHard() {super();}
    public BoundaryHard(Phase p) {super(p);}
    public BoundaryHard(Phase p, double lx, double ly) {super(p,lx,ly);}
    public void nearestImage(Space.Vector dr) {}
    public void centralImage(Space.Vector r) {}
    public void nearestImage(Vector dr) {}
    public void centralImage(Vector r) {}
    public void centralImage(Coordinate c) {}
    public void setCollisionRadius(double d) {collisionRadius = d;}
    public double getCollisionRadius() {return collisionRadius;}

    //the pair passed to this method has both atoms pointing to the same atom
    //use atom1
    public double collisionTime(AtomPair pair) {
        Atom a = pair.atom1();
        Vector r = (Vector)a.coordinate().position();
        Vector p = (Vector)a.coordinate().momentum();
        double tx = (p.x > 0.0) ? (dimensions.x - r.x - collisionRadius)/(p.x*a.rm()) : (-r.x + collisionRadius)/(p.x*a.rm());
        double ty = (p.y > 0.0) ? (dimensions.y - r.y - collisionRadius)/(p.y*a.rm()) : (-r.y + collisionRadius)/(p.y*a.rm());
        return Math.min(tx,ty);
    }
    public void bump(AtomPair pair) {
        Atom a = pair.atom1();
        Vector r = (Vector)a.coordinate().position();
        Vector p = (Vector)a.coordinate().momentum();
        double dx = (p.x > 0.0) ? Math.abs(dimensions.x - r.x - collisionRadius) : Math.abs(-r.x + collisionRadius);
        double dy = (p.y > 0.0) ? Math.abs(dimensions.y - r.y - collisionRadius) : Math.abs(-r.y + collisionRadius);
//            double dx = Math.abs(r.x/dimensions.x-0.5);   //determine which component is further from center
//            double dy = Math.abs(r.y/dimensions.y-0.5);
//            if(dx > dy) {
        if(dx < dy) {
            pAccumulator += 2*Math.abs(p.x);
            p.x *= -1;
        }
        else {
            pAccumulator += 2*Math.abs(p.y);
            p.y *= -1;
        }
    }
    //not yet implemented
    public double lastCollisionVirial() {return 0.0;}
    public Space.Tensor lastCollisionVirialTensor() {return Tensor.ZERO;}

    public double pAccumulator = 0.0;
}
