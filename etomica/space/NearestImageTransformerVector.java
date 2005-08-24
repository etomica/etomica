package etomica.space;




/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Jan 18, 2005 by kofke
 */
public class NearestImageTransformerVector implements NearestImageTransformer, java.io.Serializable {

    /* (non-Javadoc)
     * @see etomica.NearestImageTransformer#nearestImage(etomica.Space.Vector)
     */
    public void nearestImage(Vector dr) {
       if(vector == null) return;
       if(doPlus) dr.PE(vector);
       else dr.ME(vector);
    }

    public void setNearestImageVector(Vector vector) {
        this.vector = vector;
    }
    
    public void setPlus(boolean b) {
        doPlus = b;
    }
    
    private Vector vector;
    private boolean doPlus = false;
}
