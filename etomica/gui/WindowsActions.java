/**
 * WindowsActions
 *
 * The WindowsActions class is responsible for creating static action listeners to the Window drop-down menu
 * of the EtomicaMenuBar.
 *
 * @author Bryan C. Mihalick
 * 8/15/00
 */

package simulate.gui;

import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;

public class WindowsActions {
    /**
     * Static action listener for the tile event which tiles all of the internal frames
     */
    public static final ActionListener TILE = new TileAction();
    
    /**
     * Static action listener for the cascade event which cascades all of the internal frames
     */
    public static final ActionListener CASCADE = new CascadeAction();
    
    /**
     * Static action listener for the next event which selects the next internal frame in sequence
     */
    public static final ActionListener NEXT = new NextAction();
    
    /**
     * Static action listener for the drag event which only draws the outline of the component being 
     * dragged.
     */
    public static final ActionListener DRAG = new DragAction();
    
    /**
     * Static class that handles the tile event
     */
    private static class TileAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
            simulate.utility.Windows.tileWindows(((javax.swing.JDesktopPane)Etomica.DesktopFrame.desktop));
        }// end of actionPerformed
    }// end of TileAction class

    /**
     * Static class that handles the cascade event
     */
    private static class CascadeAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
            simulate.utility.Windows.cascadeWindows(((javax.swing.JDesktopPane)Etomica.DesktopFrame.desktop));
        }// end of actionPerformed
    }// end of CascadeAction class
    
    /**
     * Static class that handles the next event
     */
    private static class NextAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
            simulate.utility.Windows.selectNextWindow(((javax.swing.JDesktopPane)Etomica.DesktopFrame.desktop));
        }// end of actionPerformed
    }// end of NextAction class

    /**
     * Static class that handles the drag event
     */
    private static class DragAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
            simulate.utility.Windows.dragOutline(((javax.swing.JDesktopPane)Etomica.DesktopFrame.desktop));
        }// end of actionPerformed
    }// end of DragAction class
}// end of WindowsActions class
    