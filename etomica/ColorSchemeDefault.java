package simulate;
import java.awt.Color;
/**
 *
 * @author David Kofke
 *
 */
 
public class ColorSchemeDefault extends ColorScheme {
    
    public ColorSchemeDefault() {super();}
    public ColorSchemeDefault(Color c) {super(c);}
        
 /**
    Return without changing atom's color, which was set at initialization
  */
    
    public final void setAtomColor(Atom a) {
        return;
//        double scale = .parentSpecies.kineticEnergy()/m.parentSpecies.getNMolecules();
//        for(Atom a=m.firstAtom; a!=m.lastAtom.getNextAtom(); a=a.getNextAtom()) {
//            a.setColor(new Color((float)Math.exp(-a.kineticEnergy()/scale),(float)0.0,(float)0.0));
//        }
    }
    
}
