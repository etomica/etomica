package etomica.graphics;
import java.awt.Color;

import etomica.Simulation;

public class ConstantsGraphic extends Object {
    
    public String getVersion() {return "01.11.20";}    
    
    private ConstantsGraphic() {}   // can't instantiate class
    
  /* Colors adopted in the web textbook on molecular simulation */
    public static final Color KHAKI = new Color(153,153,102);
    public static final Color DARK_KHAKI = new Color(102,102,051);
    public static final Color BRIGHT_RED = new Color(153,000,000);
    public static final Color DARK_RED = new Color(102,000,000);
    public static final Color BLUSH = new Color(153,102,102);
    public static final Color TAN = new Color(204,204,153);
    public static final Color randomColor() {
        return new Color(Simulation.random.nextFloat(),Simulation.random.nextFloat(),Simulation.random.nextFloat());
    }
}//end of ConstantsGraphic

    
    
    