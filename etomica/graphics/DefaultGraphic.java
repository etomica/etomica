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
    
	public static Color BACKGROUND_COLOR = ConstantsGraphic.KHAKI.brighter().brighter();
	public static Color CONTRAST_COLOR = ConstantsGraphic.DARK_RED;
	public static Color BORDER_COLOR = CONTRAST_COLOR;
	public static Color SLIDER_COLOR = BACKGROUND_COLOR;//ConstantsGraphic.TAN;
	public static Color TAB_COLOR = ConstantsGraphic.DARK_RED;
	public static Color TAB_TEXT_COLOR = BACKGROUND_COLOR;
	public static Color BUTTON_COLOR = TAB_COLOR;
	public static Color BUTTON_TEXT_COLOR = TAB_TEXT_COLOR;
	public static Color PANEL_COLOR = ConstantsGraphic.TAN;
    
    
}