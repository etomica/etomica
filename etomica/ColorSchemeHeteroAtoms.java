package simulate;
import java.awt.*;

/**
 *
 * @author David Kofke
 *
 */
 
public class ColorSchemeHeteroAtoms extends ColorScheme {

    public Color[] atomColors;
    public int nAtoms;
    
    public ColorSchemeHeteroAtoms() {
        super();
    }
    public ColorSchemeHeteroAtoms(Color c) {
        super(c);
    }
    
    // need to find a way to set nAtoms to correspond to species nAtomsPerMolecule value
    public void setNAtoms(int n) {
        nAtoms = n;
        atomColors = new Color[n];
        for(int i=0; i<n; i++) {
            atomColors[i] = baseColor;
        }
    }
    public int getNAtoms() {return nAtoms;}
        
    public void setAtomColor(Atom a) {a.setColor(atomColors[a.atomIndex()]);}
    
    public void setAtomColors(int i, Color c) {atomColors[i] = c;}
    public Color getAtomColors(int i) {return atomColors[i];}
    
}
