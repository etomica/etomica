package etomica.space;

/**
 * Nearest image transformation that adds or subtracts a fixed vector amount to
 * the separation.
 */
public class NearestImageTransformerVector implements NearestImageTransformer, java.io.Serializable {

    /**
     * Default constructor leaves nearest image vector null and setPlus false.
     *
     */
    public NearestImageTransformerVector() {
    }

    /**
     * Modifies the given vector by adding (if setPlus was set true) or subtracting
     * (setPlus fals) the previously specified nearest-image vector.  If no vector
     * was previously specified, the given vector is unchanged.
     */
    public void nearestImage(Vector dr) {
        if (vector == null)
            return;
        if (doPlus)
            dr.PE(vector);
        else
            dr.ME(vector);
    }

    /**
     * Sets the vector to be added/subtracted via the nearestImage method.
     * The vector is not copied; the given instance is used, so any changes
     * to it will be reflected in the operation of nearestImage.  A null
     * vector may be specified, and will cause nearestImage to have no effect
     * when invoked.
     */
    public void setNearestImageVector(Vector vector) {
        this.vector = vector;
    }

    /**
     * Indicates if the nearest-image vector is to be added (true) or
     * substracted (false) when nearestImage is invoked.
     */
    public void setPlus(boolean b) {
        doPlus = b;
    }
    
    /**
     * Returns the current value of the "plus" flag.
     */
    public boolean isPlus() {
        return doPlus;
    }

    private Vector vector;
    private boolean doPlus = false;
    private static final long serialVersionUID = 1L;
}
