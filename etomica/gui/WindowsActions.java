/**
 * WindowsActions
 *
 * The WindowsActions class is responsible for creating static action listeners to the Window drop-down menu
 * of the EtomicaMenuBar.
 *
 * @author Bryan C. Mihalick
 * 8/15/00
 */

package etomica.gui;

import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import javax.swing.JDesktopPane;
import javax.swing.JInternalFrame;
import java.beans.PropertyVetoException;

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
     * Static class that handles the tile event and tiles all internal frames on the content pane
     */
    private static class TileAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
            JDesktopPane desktop = ((JDesktopPane)Etomica.DesktopFrame.desktop);
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
        }//end of actionPerformed method
    }// end of TileAction class

    /**
     * Static class that handles the cascade event and cascades all internal frames on the desktop
     */
    private static class CascadeAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
            JDesktopPane desktop = ((JDesktopPane)Etomica.DesktopFrame.desktop);
            JInternalFrame[] frames = desktop.getAllFrames();
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
        }// end of actionPerformed
    }// end of CascadeAction class
    
    /**
     * Static class that handles the next event and puts the focus on the next window in sequence on the desktop
     */
    private static class NextAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
            JDesktopPane desktop = ((JDesktopPane)Etomica.DesktopFrame.desktop);
            JInternalFrame[] frames = desktop.getAllFrames();
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
        }// end of actionPerformed
    }// end of NextAction class

    /**
     * Static class that handles the drag event and allows the content pane to only draw the outline 
     * of the internal frames or component being drawn on the content pane
     */
    private static class DragAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
            ((JDesktopPane)Etomica.DesktopFrame.desktop).putClientProperty("JDesktopPane.dragMode", "outline");
        }// end of actionPerformed
    }// end of DragAction class
}// end of WindowsActions class
    