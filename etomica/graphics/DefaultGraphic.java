package etomica.graphics;

import java.awt.Color;
import etomica.*;
import etomica.units.*;

public final class DefaultGraphic {
    
    static {
        Default.IS_GRAPHIC = true;
    }
        
    public static String IMAGE_DIRECTORY = "file:/" + Default.WORKING_DIRECTORY + "images/";
    
    public static String HELP_FILE = "http://www.ccr.buffalo.edu/etomica/help.html";
        
    public static String JAVADOC_FILE = "http://www.ccr.buffalo.edu/etomica/JavaDoc/index.html";

    public static /*final*/ boolean DISPLAY_USE_OPENGL = true;
    
    public static Color ATOM_COLOR = Color.black;
    
}