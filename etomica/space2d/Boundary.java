package simulate.space2D;
import simulate.Space;
import simulate.Phase;

public abstract class Boundary extends Space.Boundary {
    public static final int NONE = 0;
    public static final int PERIODIC_SQUARE = 1;
    public static final int HARD = 2;
    public static final int SLIDING_BRICK = 3;
    public static final String[] TAGS = {"None", "Periodic Square", "Hard", "Sliding Brick"};
    public Boundary() {super();}
    public Boundary(Phase p) {super(p);}
    public abstract void nearestImage(Vector dr);
    public abstract void centralImage(Vector r);
    public abstract void centralImage(Coordinate c);
}

