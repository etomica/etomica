package simulate;
import java.awt.*;

/**
 *
 * @author David Kofke
 *
 */
 
public abstract class ColorScheme extends Component {

    public Color baseColor;
    
    public ColorScheme() {
        this(Color.black);
    }
    public ColorScheme(Color c) {
        setBaseColor(c);
    }
        
    public abstract void setAtomColor(Atom a);
    
    public void initializeAtomColor(Atom a) {
        a.setColor(baseColor);
    }
    
    public void initializeMoleculeColor(Molecule m) {
        for(Atom a=m.firstAtom(); a!=m.terminationAtom(); a=a.getNextAtom()) {
            initializeAtomColor(a);
        }
    }
    
    public final void setBaseColor(Color c) {baseColor = c;}
    public final Color getBaseColor() {return baseColor;}
}
