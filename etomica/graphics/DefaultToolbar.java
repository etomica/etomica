package etomica.graphics;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;


public class DefaultToolbar {

	private JMenuBar mBar;
	private String[] DEFAULT_AUTHORS = {"Dr. David A. Kofke", "Dr. Andrew Schultz" };
    private String[] DEFAULT_CONTRIB = { "Robert Rassler" };
    private String[] authors = DEFAULT_AUTHORS;
    private String[] contributors = DEFAULT_CONTRIB;

    public DefaultToolbar() {
    	this(null, "");
    }

    public DefaultToolbar(final JPanel parent, String appName) {
    	mBar = new JMenuBar();

    	
    	// File menu
    	JMenu fileMenu = new JMenu("File");

    	JMenuItem exitBtn = new JMenuItem("Exit");
    	exitBtn.addActionListener(new ActionListener() {
    		public void actionPerformed(ActionEvent ev) {
    			System.exit(0);
    		}
    	});
    	fileMenu.add(exitBtn);

    	// Help menu
    	JMenu helpMenu = new JMenu("Help");
    	final String aboutString = "About " + appName;
    	JMenuItem aboutBtn = new JMenuItem(aboutString);
    	aboutBtn.addActionListener(new ActionListener() {
    		public void actionPerformed(ActionEvent ev) {
    			AboutBoxWindow about =
    				new AboutBoxWindow(parent,
    					               aboutString,
    					               authors,
    					               contributors);
    			about.setVisible(true);
    		}
    	});

    	aboutBtn.setEnabled(true);
    	helpMenu.add(aboutBtn);

    	mBar.add(fileMenu);
    	mBar.add(helpMenu);

    }

    public JMenuBar graphic() {
    	return mBar;
    }

}
