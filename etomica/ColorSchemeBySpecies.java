//includes a main method to test and demonstrate class
package etomica;
import java.awt.Color;

/**
* Color scheme that sets and uses a different color scheme for each species.
*
* @author David Kofke
*/
public class ColorSchemeBySpecies extends ColorScheme {
    
    public static String getVersion() {return "01.01.17/"+ColorScheme.getVersion();}
    private SpeciesColor[] speciesColor = new SpeciesColor[0];
    
    public ColorSchemeBySpecies() {super();}
    
    /**
     * Add species to list of species affected by this color scheme
     * and defines particular scheme to be used to color the species.
     */
     //what if speciesIndex changes after adding to this?
    public void addSpecies(Species s, ColorScheme cs) {
        if(s == null || cs == null) return;
        if(s.index() >= speciesColor.length) {//make a bigger array
            SpeciesColor[] temp = new SpeciesColor[s.index()+1];
            for(int i=0; i<speciesColor.length; i++) {
                temp[i] = speciesColor[i];
            }
            speciesColor = temp;
        }
        if(getPhase() != null) cs.setIterator(s.makeAtomIterator(getPhase()));
        speciesColor[s.index()] = new SpeciesColor(s, cs);
    }
    
    /**
     * Override of superclass method to make species iterators for all added species
     */
    public void setPhase(Phase p) {
        super.setPhase(p);
        for(int i=0; i<speciesColor.length; i++) {
            speciesColor[i].colorScheme.setIterator(
                speciesColor[i].species.makeAtomIterator(getPhase()));
        }
    }

    /**
     * Colors the atoms according to the scheme for its species
     */
    public void colorAtom(Atom a) {  
        speciesColor[a.parentSpecies().index()].colorScheme.colorAtom(a);
    }
    
    /**
     * Colors all atoms by looping over all species, invoking the color scheme for each
     */
    public void colorAllAtoms() {
        for(int i=0; i<speciesColor.length; i++) {
            speciesColor[i].colorScheme.colorAllAtoms();
        }
    }

    /**
     * Inner class used to join a species and a color scheme
     */
    private static final class SpeciesColor {
        Species species;
        ColorScheme colorScheme;
        SpeciesColor(Species s, ColorScheme c) {
            species = s;
            colorScheme = c;
        }
    }
    
    /**
     * Method to test and demonstrate use of this class
     */
    public static void main(String[] args) {

        Simulation.instance = new Simulation(new Space2D());
        
        //set up a system of three species
	    IntegratorHard integratorHard1 = new IntegratorHard();
	    Controller controller1 = new Controller();
	    Phase phase1 = new Phase();
	    SpeciesDisks speciesDisks0 = new SpeciesDisks(10,3); //nM, nA
	    SpeciesDisks speciesDisks1 = new SpeciesDisks(5,1);
	    SpeciesDisks speciesDisks2 = new SpeciesDisks(5,1);
	    P2Tether = new P2Tether();
	    P2HardDisk p2HardDisk00 = new P2HardSphere();  
	    P2HardDisk p2HardDisk01 = new P2HardSphere();  
	    P2HardDisk p2HardDisk02 = new P2HardSphere();  
	    P2HardDisk p2HardDisk11 = new P2HardSphere(); 
	    P2HardDisk p2HardDisk12 = new P2HardSphere();  
	    P2HardDisk p2HardDisk22 = new P2HardSphere();
	    Simulation.instance.elementCoordinator().go();
        Potential2.Agent potentialAgent = (Potential2.Agent)p2HardDisk00.getAgent(phase);
        potentialAgent.setIterator(new AtomPairIterator(phase1, speciesDisks0, speciesDisks0));
	    
	  /*  P2HardDisk01.setSpeciesIndex(0,1);
	    P2HardDisk02.setSpeciesIndex(0,2);
	    P2HardDisk11.setSpeciesIndex(1,1);
	    P2HardDisk12.setSpeciesIndex(1,2);
	    P2HardDisk22.setSpeciesIndex(2,2);
*/
	    DisplayPhase displayPhase1 = new DisplayPhase();
	    IntegratorMD.Timer timer = integratorHard1.new Timer(integratorHard1.new ChronoMeter());
	    timer.setUpdateInterval(10);
		Simulation.instance.setBackground(Color.yellow);
		
		//set up color scheme
		ColorSchemeBySpecies colorScheme = new ColorSchemeBySpecies();
//	for debugging	ColorSchemeHeteroAtoms cs0 = new ColorSchemeHeteroAtoms();
		ColorSchemeHeteroAtoms cs0 = new ColorSchemeNull();
		ColorSchemeTemperature cs1 = new ColorSchemeTemperature();
		ColorSchemeNull cs2 = new ColorSchemeNull();
		colorScheme.addSpecies(speciesDisks0, cs0);
		colorScheme.addSpecies(speciesDisks1, cs1);
		colorScheme.addSpecies(speciesDisks2, cs2);
		displayPhase1.setColorScheme(colorScheme);
		cs0.setNAtoms(speciesDisks0.getAtomsPerMolecule());
        cs0.setAtomColors(new Color[] {Color.yellow, Color.green, Color.cyan});
		
		
		Simulation.instance.elementCoordinator.go(); 
		Simulation.makeAndDisplayFrame(Simulation.instance);                                    
    }
        
} //end of ColorSchemeBySpecies
