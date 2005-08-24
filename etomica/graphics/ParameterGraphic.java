package etomica.graphics;
import etomica.atom.Parameter;

public abstract class ParameterGraphic extends Parameter {
    
    public interface Color {
        public java.awt.Color getColor();
        public void setColor(java.awt.Color c);
    }
}
    