package etomica.graphics;
import java.awt.Color;

public class ParameterColor extends ParameterGraphic implements ParameterGraphic.Color {
    
    private java.awt.Color color = DefaultGraphic.ATOM_COLOR;
    
    public java.awt.Color getColor() {return color;}
    public void setColor(java.awt.Color c) {color = c;}
    
}
    