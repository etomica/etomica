package etomica.graphics;
import etomica.*;

/**
 * Changes the color of an atom to some specified value.
 *
 * Not implemented.
 */
public class AtomActionChangeColor extends AtomAction implements Action.Undoable {
    
    private Color newColor = Color.blue;
    private Color oldColor = Color.black;
    public AtomActionChangeColor() {this(Color.red);}
    public AtomActionChangeColor(Color c) {setNewColor(c);}
    public void setNewColor(Color c) {newColor = c;}
    public Color getNewColor() {return newColor;}
    public void actionPerformed(Atom a) {
        atom = a;
 /*       oldColor = a.getColor();
        a.setColor(newColor);*/
    }
    public void attempt() {
        if(atom != null) atom.setColor(newColor);
    }
    public void undo() {
        if(atom != null) atom.setColor(oldColor);
    }
}//end of AtomActionChangeColor
    
