//includes a main method to demonstrate use and to test
package etomica;
import java.awt.Color;

/**
 * Color scheme that permits different colors for each atom of a molecule.
 *
 * @author David Kofke
 */
 
 //needs work to tie better with species
public class ColorSchemeHeteroAtoms extends ColorScheme {

    private Color[] atomColors;
    private int nAtoms;
    private Species species;
    
    public ColorSchemeHeteroAtoms() {
        super();
    }
    public ColorSchemeHeteroAtoms(Color c) {
        super(c);
    }
    
    public void setSpecies(Species s) {
        species = s;
        setNAtoms(species.getAtomsPerMolecule());
    }
    public Species getSpecies() {return species;}
    
    public void setNAtoms(int n) {
        nAtoms = n;
        atomColors = new Color[n];
        for(int i=0; i<n; i++) {
            atomColors[i] = baseColor;
        }
    }
    public int getNAtoms() {return nAtoms;}
     
        //does not check for array out of bounds; need a stronger tie to the species
    public void colorAtom(Atom a) {a.setColor(atomColors[a.atomIndex()]);}
    
//    public void colorAllAtoms() {}
    
    public void setAtomColors(Color[] colors) {
        if(colors.length != nAtoms) {
            System.out.println("Error in ColorSchemeHeteroAtoms.setAtomColors");
            System.exit(1);
        }
        atomColors = colors;
    }
    public Color[] getAtomColors() {return atomColors;}
    public void setAtomColors(int i, Color c) {
        atomColors[i] = c;
        initialColorAllAtoms();
    }
    public Color getAtomColors(int i) {return atomColors[i];}
    
    /**
     * Demonstrates how this class is implemented.  
     * Appropriate for system having only one species.  For multiple species, see ColorSchemeBySpecies
     */
    public static void main(String[] args) {
        etomica.simulations.HSMD2D sim = new etomica.simulations.HSMD2D();
        Simulation.instance = sim;

        //part unique to this example
             //get handles to components we need
        SpeciesDisks species = sim.species;
        DisplayPhase display = sim.display;
             //set species to have 3 atoms per molecule to use this color scheme
        species.setAtomsPerMolecule(3);    
        P1TetherHardDisk potentialTether = new P1TetherHardDisk(); //an intramolecular potential
             //instantiate color scheme and link appropriately
        ColorSchemeHeteroAtoms colorScheme = new ColorSchemeHeteroAtoms();
        colorScheme.setSpecies(species); //sets the number of atoms
        colorScheme.setAtomColors(new Color[] {Color.red, Color.yellow, Color.blue});
        display.setColorScheme(colorScheme);
        //end of unique part
		sim.elementCoordinator.go(); 
        Simulation.makeAndDisplayFrame(sim);         
    }//end main
}
