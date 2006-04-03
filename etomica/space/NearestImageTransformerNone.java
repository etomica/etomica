package etomica.space;

/**
 * Implementation of NearestImageTransformer that does nothing.  Suitable when
 * periodic boundaries are not involved.
 * @author Andrew Schultz
 */
public class NearestImageTransformerNone implements NearestImageTransformer {

    /**
     * Does nothing
     */
    public void nearestImage(Vector dr) {
    }

}
