/**
 * ButtonActions
 *
 * The ButtonActions class is responsible for creating static action listeners to the EtomicaToolBar
 * buttons.
 *
 * @author Bryan C. Mihalick
 * 8/14/00
 */
 
package etomica.gui;

import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import javax.swing.JFileChooser;
import javax.swing.JEditorPane;
import javax.swing.JScrollPane;
import javax.swing.JInternalFrame;
import java.awt.Component;
import com.symantec.itools.javax.swing.JToolBarSeparator;
import com.symantec.itools.javax.swing.icons.ImageIcon;

public class ButtonActions {
    
    public static final String version() {return "01.03.04.0";}

    /**
     * Static action listener that pastes a copied component to the FormDesign
     */
    public static final ActionListener PASTE = new PasteAction();
    
    /**
     * Static action listener that cuts a component from the FormDesign
     */
    public static final ActionListener CUT = new CutAction();
    
    /**
     * Static action listener that copies a component from the FormDesign
     */
    public static final ActionListener COPY = new CopyAction(); 
    
    /**
     * This will eventually call a paste method for pasting a copy of a simulation component
     */
    private static class PasteAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {

        }// end of actionPerformed
    }// end of PasteAction
    
    /**
     * This will eventually call a cut method for deleting an instance of a simulation component
     */
    private static class CutAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
            
        }// end of actionPerformed
    }// end of CutAction

    /**
     * This will eventually call a copy method for saving a copy of an instance of a simulation component
     */
    private static class CopyAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {

        }// end of actionPerformed
    }// end of CopyAction    
}// end of ButtonActions class
    