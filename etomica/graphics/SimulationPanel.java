package etomica.graphics;
import javax.swing.JPanel;

/**
 * Main panel used to display the graphical elements of a simulation.  
 * An instance of this panel is first constructed when the panel() method
 * of a Simulation instance is invoked.  This method is called when the panel
 * is added to an applet or a window, and it is also invoked by the default
 * Mediator class when a Display or Device are added to a Simulation.
 *
 * @author David Kofke
 */

/* History
 * 
 * 01/03/03 (DAK) added constructor that omits addition of internal panels
 */
 
public class SimulationPanel extends JPanel {
    
	public final javax.swing.JTabbedPane displayPanel = new javax.swing.JTabbedPane();
	public final javax.swing.JPanel displayBoxPanel = new JPanel(new java.awt.GridBagLayout());
//    public final javax.swing.JPanel devicePanel = new JPanel(new java.awt.GridLayout(0,1),false);
    public final javax.swing.JPanel devicePanel = new JPanel(new java.awt.GridBagLayout());
//    public final javax.swing.JPanel devicePanel = new JPanel(new java.awt.FlowLayout());

    
    public SimulationPanel() {
        add(devicePanel);
        add(displayBoxPanel);
        add(displayPanel);
        setSize(800,550);
        setLayout(new java.awt.FlowLayout());
//        setBackground(DefaultGraphic.BACKGROUND_COLOR);
/*        setLayout(new java.awt.BorderLayout());
        add(devicePanel, java.awt.BorderLayout.NORTH);
        add(displayPanel, java.awt.BorderLayout.EAST);
        add(displayBoxPanel, java.awt.BorderLayout.WEST);*/
/*        setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gbc = new java.awt.GridBagConstraints();
        gbc.gridx = 0;
        add(devicePanel, gbc);
        add(displayBoxPanel, gbc);
        gbc.gridx = 1;
        add(displayPanel, gbc);*/
        //workaround for JTabbedPane bug in JDK 1.2
        displayPanel.addChangeListener(
            new javax.swing.event.ChangeListener() {
                public void stateChanged(javax.swing.event.ChangeEvent event) {
                    displayPanel.invalidate();
                    displayPanel.validate();
                }
        });
        
 /*       displayPanel.addMouseListener( 
            new java.awt.event.MouseAdapter() {
                public void mouseClicked(java.awt.event.MouseEvent evt) {
                    System.out.println("Click");
                }
            });*/
    }
}//end of SimulationPanel