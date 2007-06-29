package etomica.box;

/**
 * Event that conveys some happening with respect to a box or the things it contains.
 *
 * @see BoxListener
 */
public class BoxEvent implements java.io.Serializable {
    
    public BoxEvent(Box box) {
        super();
        this.box = box;
    }
    
    public Box getBox() {
        return box;
    }
    
    private final Box box;
    private static final long serialVersionUID = 1L;
}
    