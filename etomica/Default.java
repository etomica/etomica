package simulate;

import java.awt.Color;
import simulate.units.*;

public final class Default {
    
    public static String CLASS_DIRECTORY = "D:/kofke/development/lib/simulate";
    
    public static String IMAGE_DIRECTORY = "D:/kofke/development/images/";
        
    public static double ATOM_SIZE = 3.0;  //Angstroms
    
    public static double ATOM_MASS = 40.0; //Daltons
    
    public static Color ATOM_COLOR = Color.black;
    
    public static int MOLECULE_COUNT = 20;
    
    public static double BOX_SIZE = 30.0;  //Angstroms
    
    public static double TEMPERATURE = Kelvin.UNIT.toSim(300.);
    
    public static double PRESSURE = 1.0;
    
    public static double POTENTIAL_WELL = Kelvin.UNIT.toSim(300.);
    
    public static double POTENTIAL_CUTOFF = 2.5; //dimensionless multiplier for cutoff distance of potential
    
    public static double TIME_STEP = 0.05;  //picoseconds
    
}