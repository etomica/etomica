package etomica.graphics;
import java.awt.Color;

import etomica.Parameter;
import etomica.atom.Atom;
import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomType;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;

/**
 * Colors the atom according to the color given by its type field.  Instantiation
 * of this class leads to the placement of a ParameterColor object with each
 * existing and subsequent instance of AtomType.  This parameter is used to store
 * and reference the color associated with each type.
 *
 * @author David Kofke
 */

public final class ColorSchemeByType extends ColorScheme implements Parameter.Source {
    
    public static final int colorIndex = AtomType.requestParameterIndex(
        new Parameter.Source() {
            public Parameter makeParameter() {return new ParameterColor();}
        });
        
    public ColorSchemeByType() {}
            
    /**
     * Initialize atom color to the color of its type
     */
    public final Color atomColor(Atom a) {
        return(((ParameterColor)a.type.parameter[colorIndex]).getColor());}
  
    //implementation of Parameter.Source interface
    public Parameter makeParameter() {return new ParameterColor();}
    
    
    public static void setColor(Species s, Color c) {
        ((ParameterColor)(s.moleculeFactory()).getType().parameter[ColorSchemeByType.colorIndex]).setColor(c);
    }
   
    public static void setColor(AtomType type, Color c) {
        ((ParameterColor)type.parameter[ColorSchemeByType.colorIndex]).setColor(c);
    }
    
    public static Color getColor(SpeciesSpheresMono s) {
        return ((ParameterColor)((AtomFactoryMono)s.moleculeFactory()).getType().parameter[ColorSchemeByType.colorIndex]).getColor();
    }
        /*
    public static void setColor(SpeciesWater s, Color c) {
        ((ParameterColor)((AtomFactoryHetero)s.moleculeFactory()).childFactory.type().parameter[ColorSchemeByType.colorIndex]).setColor(c);
    }
    */       
   
    /**
     * Demonstrates how this class is implemented.
     */
//    public static void main(String[] args) {
//        Simulation.instance = new SimulationGraphic(new Space2D());
//	    IntegratorHard integratorHard = new IntegratorHard();
//	    SpeciesSpheresMono speciesBlue = new SpeciesSpheresMono();
//	    SpeciesSpheresMono speciesRed = new SpeciesSpheresMono();
//	    
//	    Phase phase = new Phase();
//	    Potential2 potential = new P2HardSphere();
//	    Controller controller = new Controller();
//	    DisplayPhase displayPhase = new DisplayPhase();
//	    
//	    //this is the special part
//        displayPhase.setColorScheme(new ColorSchemeByType());        
//        ((ParameterColor)((AtomFactoryMono)speciesBlue.moleculeFactory()).type().parameter[ColorSchemeByType.colorIndex]).setColor(Color.blue);
//        ((ParameterColor)((AtomFactoryMono)speciesRed.moleculeFactory()).type().parameter[ColorSchemeByType.colorIndex]).setColor(Color.red);
//        //--------------------
//                
//        //this method call invokes the mediator to tie together all the assembled components.
//		Simulation.instance.elementCoordinator.go();
//		                                    
//		((SimulationGraphic)Simulation.instance).panel().setBackground(java.awt.Color.yellow);
//        SimulationGraphic.makeAndDisplayFrame(Simulation.instance);
//    }//end of main
    
}
