package etomica.graphics;
import etomica.*;

public interface Drawable {
    public void draw(java.awt.Graphics g, int[] origin, double scale);
}
