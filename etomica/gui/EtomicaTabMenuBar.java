/**
 * EtomicaTabMenuBar
 *
 * This class creates the tab menu bar for the Etomica simulation environment.  This includes, creating
 * the JTabbedPane itself, giving the pane and its tabs titles, adding the tabs, and adding components 
 * to the tabs.
 *
 * @author Bryan C. Mihalick
 * 8/14/00
 */

package etomica.gui;

import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JComponent;
import java.awt.Image;
import java.awt.Toolkit;

public class EtomicaTabMenuBar extends javax.swing.JToolBar {
    /**
     * Main pane that contains all tabs
     */
    javax.swing.JTabbedPane tabMenu = new javax.swing.JTabbedPane(javax.swing.SwingConstants.BOTTOM);

    /**
     * Images displayed on buttons so that they can be distinguished from one another
     */
    javax.swing.ImageIcon jButton = new javax.swing.ImageIcon("javax.swing.JButton.ico");
    javax.swing.ImageIcon jFrame = new javax.swing.ImageIcon("javax.swing.JFrame.ico");
    javax.swing.ImageIcon jPanel = new javax.swing.ImageIcon("javax.swing.JPanel.ico");
    javax.swing.ImageIcon jSlider = new javax.swing.ImageIcon("javax.swing.JSlider.ico");
    javax.swing.ImageIcon jTextField = new javax.swing.ImageIcon("javax.swing.JTextField.ico");

    /**
     * Buttons added to each tab of the JTabbedPane
     */
    javax.swing.JButton button1 = new javax.swing.JButton("JButton");
    javax.swing.JButton button2 = new javax.swing.JButton("JFrame");
    javax.swing.JButton button3 = new javax.swing.JButton("JPanel");
    javax.swing.JButton button4 = new javax.swing.JButton("JSlider");
    javax.swing.JButton button5 = new javax.swing.JButton("JTextField");

    /**
     * Images displayed on buttons so that they can be distinguished from one another
     */
    Image image1 = Toolkit.getDefaultToolkit().getImage("javax.swing.JButton.ico");
    Image image2 = Toolkit.getDefaultToolkit().getImage("javax.swing.JFrame.ico");
    Image image3= Toolkit.getDefaultToolkit().getImage("javax.swing.JPanel.ico");
    Image image4 = Toolkit.getDefaultToolkit().getImage("javax.swing.JSlider.ico");
    Image image5 = Toolkit.getDefaultToolkit().getImage("javax.swing.JTextField.ico");
    
    /**
     * Constructor that makes the tab menu bar for the Etomica environment.  A tab is made for each type
     * of simulation component.  Each of these tabs contains buttons that allow an instance of each
     * simulation component type to be added to the FormDesign
     */
    public EtomicaTabMenuBar() {
        this.setBounds(200, 0, 575, 56);

        /**
         * Add tabs with the superclass name of the respective simulation components that they contain 
	     */
	    tabMenu.addTab("Phase", new javax.swing.JPanel());
	    tabMenu.addTab("Potentials", new javax.swing.JPanel());
	    tabMenu.addTab("Species", new javax.swing.JPanel());
	    tabMenu.addTab("Integrators", new javax.swing.JPanel());
	    tabMenu.addTab("Controllers", new javax.swing.JPanel());
	    tabMenu.addTab("Displays", new javax.swing.JPanel());
        tabMenu.addTab("Meters", new javax.swing.JPanel());
	    tabMenu.addTab("Devices", new javax.swing.JPanel());
	    tabMenu.addTab("Actions", new javax.swing.JPanel());
	    
        /**
	     * Create dummy button icons to test displaying ability of the JTabbedPane tabs
	     */
	    button1.setIcon(jButton);
	    button2.setIcon(jFrame);
	    button3.setIcon(jPanel);
	    button4.setIcon(jSlider);
	    button5.setIcon(jTextField);
	    
        /**
	     * Add the dummy button components to the tabs with their respective icons
	     */
	    for(int i = 0; i < tabMenu.getComponentCount(); i++){
	        ((JComponent)tabMenu.getComponent(i)).add(new JButton(new ImageIcon("icons.Simulation_COLOR_16x16"))); 
	        ((JComponent)tabMenu.getComponent(i)).add(new JButton(new ImageIcon(image2)));
	        ((JComponent)tabMenu.getComponent(i)).add(new JButton(new ImageIcon(image3)));
	        ((JComponent)tabMenu.getComponent(i)).add(new JButton(new ImageIcon(image4)));
	        ((JComponent)tabMenu.getComponent(i)).add(new JButton(new ImageIcon(image5)));
	    }
	    
	    this.add(tabMenu);
	}// end of EtomicaTabMenuBar constructor
}// end of EtomicaTabMenuBar class