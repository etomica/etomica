/**
 * Windows
 *
 * The Windows class is responsible for creating static methods to the Window drop-down menu action
 * listeners of the EtomicaMenuBar (these are located in WindowsActions).
 *
 * @author Bryan C. Mihalick
 * 8/15/00
 */

package etomica.utility;

import javax.swing.JDesktopPane;
import javax.swing.JInternalFrame;
import java.beans.PropertyVetoException;

public class Windows {
    /**
     * Static method that tiles all internal frames on the content pane
     */
    public static void tileWindows(JDesktopPane desktop){
        JInternalFrame[] frames = desktop.getAllFrames();

        // count frames that aren't iconized
        int frameCount = 0;
        for (int i = 0; i < frames.length; i++){
            if (!frames[i].isIcon())
                frameCount++;
        }

        int rows = (int)Math.sqrt(frameCount);
        int cols = frameCount / rows;
        int extra = frameCount % rows;
            // number of columns with an extra row

        int width = desktop.getWidth() / cols;
        int height = (desktop.getHeight()-desktop.getComponentAt(0,0).getHeight()) / rows;
        int r = 0;
        int c = 0;
        for (int i = 0; i < frames.length; i++) {
            if (!frames[i].isIcon()) {
                try {
                    frames[i].setMaximum(false);
                    frames[i].reshape(c * width, r * height + 
                        desktop.getComponentAt(0,0).getHeight(), width, height);
                    r++;
                    if (r == rows) {
                        r = 0;
                        c++;
                        if (c == cols - extra) {
                            // start adding an extra row
                            rows++;
                            height = (desktop.getHeight()-desktop.getComponentAt(0,0).getHeight()) / rows;
                        }
                    }
                }
                catch(PropertyVetoException e){}
            }
        }
    }//end of tileWindows method
    
    /**
     * Static method that cascades all internal frames on the content pane
     */
    public static void cascadeWindows(JDesktopPane desktop)
    {  JInternalFrame[] frames = desktop.getAllFrames();
        int x = 0;
        int y = desktop.getComponentAt(0,0).getHeight();
        int width = desktop.getWidth() / 2;
        int height = (desktop.getHeight()-desktop.getComponentAt(0,0).getHeight()) / 2;

        for (int i = 0; i < frames.length; i++) {
            if (!frames[frames.length-i-1].isIcon()) {
                try {
                    /* try to make maximized frames resizable
                     *  this might be vetoed
                     */
                    frames[frames.length-i-1].setMaximum(false);
                    frames[frames.length-i-1].reshape(x, y, width, height);

                    x += 15;
                    y += 15;
                    
                    // wrap around at the desktop edge
                    if (x + width > desktop.getWidth()) x = 0;
                    if (y + height > (desktop.getHeight()-desktop.getComponentAt(0,0).getHeight())) y = desktop.getComponentAt(0,0).getHeight();
                }
                catch(PropertyVetoException e){}
            }
        }
    }//end of cascadeWindows method

    /**
     * Static method that selects the next window in sequence on the content pane
     */
    public static void selectNextWindow(JDesktopPane desktop)
    {  JInternalFrame[] frames = desktop.getAllFrames();
        for (int i = 0; i < frames.length; i++) {
            if (frames[i].isSelected()){
               /* find next frame that isn't an icon and can be
                selected
                */
                try {
                    int next = i + 1;
                    while (next != i && frames[next].isIcon())
                        next++;
                    if (next == i) return;
                    // all other frames are icons or veto selection
                    frames[next].setSelected(true);
                    frames[next].toFront();
                    return;
                }
                catch(PropertyVetoException e){}
            }
        }
    }// end of selectNextWindow method

    /**
     * Static method that allows the content pane to only draw the outline of the internal frames or 
     * component being drawn on the content pane
     */
    public static void dragOutline(JDesktopPane desktop) {
        
        desktop.putClientProperty("JDesktopPane.dragMode", "outline");
    }// end of dragOutline
}// end of Windows class
    
    
    

