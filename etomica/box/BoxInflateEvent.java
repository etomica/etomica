package etomica.box;


/**
 * Event that conveys that the boundary dimensions of a box have changed.
 */
public class BoxInflateEvent extends BoxEvent {

    public BoxInflateEvent(Box box) {
        super(box);
    }

    private static final long serialVersionUID = 1L;
}
