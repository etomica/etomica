/**
 * HelpActions
 *
 * The HelpActions class is responsible for creating static action listeners to the Help drop-down menu
 * of the EtomicaMenuBar.
 *
 * @author Bryan C. Mihalick
 * 8/14/00
 */

package simulate.gui;

import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.awt.Component;
import javax.swing.JEditorPane;
import javax.swing.JInternalFrame;
import javax.swing.JScrollPane;
import java.io.IOException;

public class HelpActions {
    /**
     * Static action listener that displays the documentation html pages for the simulation classes
     */
    public static final ActionListener HELP = new HelpAction();
    
    /**
     * Static action listener that displays information about the Etomica environment
     */
    public static final ActionListener ABOUT = new AboutAction();
    
    /**
     * Handles the Help event and displays the simulation documentation in a web browser internal frame
     */
    private static class HelpAction implements ActionListener {
        private int nextFrameX;
        private int nextFrameY;
        private int frameDistance;
        
        public void actionPerformed(ActionEvent event) {
            try{  
                java.net.URL fileUrl = new java.net.URL("file:\\winnt\\profiles\\mihalick\\personal\\molsim\\semest~1\\api\\packages.html");
                createInternalFrame(createEditorPane(fileUrl),"Packages");
            }
            catch(java.net.MalformedURLException e){}
        }// end of actionPerformed
        
        /**
         * Create an editor pane that follows hyperlink clicks
         */
        public Component createEditorPane(java.net.URL u){  
            JEditorPane editorPane = new JEditorPane();
            editorPane.setEditable(false);
            editorPane.addHyperlinkListener(new javax.swing.event.HyperlinkListener()
                {  public void hyperlinkUpdate(javax.swing.event.HyperlinkEvent event)
                    {   if (event.getEventType() == javax.swing.event.HyperlinkEvent.EventType.ACTIVATED){
                            createInternalFrame(createEditorPane(event.getURL()),
                            event.getURL().toString());
                        }
                    }
                });
            try
            {  editorPane.setPage(u);
            }
            catch(java.io.IOException e)
            {  editorPane.setText("Error: " + e);
            }
            return new JScrollPane(editorPane);
        }// end of createEditorPane
        
        /**
         * Create an internal frame to contain the editor pane
         */
        public void createInternalFrame(Component c, String t){  
            JInternalFrame iframe = new JInternalFrame(t,
                true, // resizable
                true, // closable
                true, // maximizable
                true); // iconifiable

            iframe.getContentPane().add(c);
            Etomica.DesktopFrame.desktop.add(iframe);

    //        iframe.setFrameIcon(new ImageIcon("document.gif"));

            // add listener to confirm frame closing
            iframe.addVetoableChangeListener(Etomica.DesktopFrame.etomicaFrame);

            // position frame
            int width = Etomica.DesktopFrame.desktop.getWidth() / 2;
            int height = Etomica.DesktopFrame.desktop.getHeight() / 2;
            iframe.reshape(nextFrameX, nextFrameY, width, height);

            iframe.show();

            // select the frame--might be vetoed
            try {  
                iframe.setSelected(true);
            }
            catch(java.beans.PropertyVetoException e){}

            /** 
             * If this is the first time, compute distance between
             * cascaded frames
             */

            if (frameDistance == 0)
                frameDistance = iframe.getHeight() - iframe.getContentPane().getHeight();

            // compute placement for next frame

            nextFrameX += frameDistance;
            nextFrameY += frameDistance;
            if (nextFrameX + width > Etomica.DesktopFrame.desktop.getWidth())
                nextFrameX = 0;
            if (nextFrameY + height > Etomica.DesktopFrame.desktop.getHeight())
                nextFrameY = 0;
        }// end of createInternalFrame
    }// end of HelpAction

    /**
     * Handles the about event
     */
    private static class AboutAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
		    try {
			    // JAboutDialog Create with owner and show as modal
				    JAboutDialog JAboutDialog1 = new JAboutDialog(Etomica.DesktopFrame.etomicaFrame);
				    JAboutDialog1.setResizable(false);
				    JAboutDialog1.setLocation(337, 284);
				    JAboutDialog1.setModal(true);
				    JAboutDialog1.show();
		    } catch (Exception e) {}
            
        }// end of actionPerformed
    }// end of AboutAction
}// end of HelpActions class
    