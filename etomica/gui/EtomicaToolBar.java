/**
 * EtomicaToolBar
 *
 * This class creates the editing tool bar for the Etomica simulation environment.  This includes, giving
 * icons, adding JButton instances, adding tooltips, and adding actionListeners that link to
 * static action listeners defined in the <JMenu name> + Actions classes.
 *
 * @author Bryan C. Mihalick
 * 8/14/00
 */

package simulate.gui;

import javax.swing.JButton;
import java.net.URL;

public class EtomicaToolBar extends javax.swing.JToolBar {
    
    /**
     * Icons displayed on editing buttons to distinguish them
     */
    com.symantec.itools.javax.swing.icons.ImageIcon newIcon = new com.symantec.itools.javax.swing.icons.ImageIcon();
	com.symantec.itools.javax.swing.icons.ImageIcon openIcon = new com.symantec.itools.javax.swing.icons.ImageIcon();
	com.symantec.itools.javax.swing.icons.ImageIcon saveIcon = new com.symantec.itools.javax.swing.icons.ImageIcon();
	com.symantec.itools.javax.swing.icons.ImageIcon cutIcon = new com.symantec.itools.javax.swing.icons.ImageIcon();
	com.symantec.itools.javax.swing.icons.ImageIcon copyIcon = new com.symantec.itools.javax.swing.icons.ImageIcon();
	com.symantec.itools.javax.swing.icons.ImageIcon pasteIcon = new com.symantec.itools.javax.swing.icons.ImageIcon();
	com.symantec.itools.javax.swing.icons.ImageIcon aboutIcon = new com.symantec.itools.javax.swing.icons.ImageIcon();

    /**
     * Editing buttons used on the toolbar along with separators to break them up
     */
	JButton newButton = new JButton();
	JButton openButton = new JButton();
	JButton saveButton = new JButton();
	com.symantec.itools.javax.swing.JToolBarSeparator JToolBarSeparator1 = new com.symantec.itools.javax.swing.JToolBarSeparator();
	JButton cutButton = new JButton();
	JButton copyButton = new JButton();
	JButton pasteButton = new JButton();
	com.symantec.itools.javax.swing.JToolBarSeparator JToolBarSeparator2 = new com.symantec.itools.javax.swing.JToolBarSeparator();
	JButton aboutButton = new JButton();
    
    /**
     * Constructor for the tool bar of the Etomica environment.  It creates all the editing buttons, adds
     * them to the toolbar, gives them icons and tooltips, and finally adds actionlisteners.
     */
    public EtomicaToolBar() {
        
        /**
         * Find and add icons to each button instance
         */
        try {
			//newIcon.setImageLocation(symantec.itools.net.RelativeURL.getURL("images/new.gif"));
			newIcon.setImageLocation(new URL(simulate.Default.IMAGE_DIRECTORY+"new.gif"));
		}
		catch (java.net.MalformedURLException error) { }
		//$$ newIcon.move(144,312);
		try {
			openIcon.setImageLocation(new URL(simulate.Default.IMAGE_DIRECTORY+"open.gif"));
		}
		catch (java.net.MalformedURLException error) { }
		//$$ openIcon.move(120,312);
		try {
			saveIcon.setImageLocation(new URL(simulate.Default.IMAGE_DIRECTORY+"save.gif"));
		}
		catch (java.net.MalformedURLException error) { }
		//$$ saveIcon.move(96,312);
		try {
			cutIcon.setImageLocation(new URL(simulate.Default.IMAGE_DIRECTORY+"cut.gif"));
		}
		catch (java.net.MalformedURLException error) { }
		//$$ cutIcon.move(72,312);
		try {
			copyIcon.setImageLocation(new URL(simulate.Default.IMAGE_DIRECTORY+"copy.gif"));
		}
		catch (java.net.MalformedURLException error) { }
		//$$ copyIcon.move(48,312);
		try {
			pasteIcon.setImageLocation(new URL(simulate.Default.IMAGE_DIRECTORY+"paste.gif"));
		}
		catch (java.net.MalformedURLException error) { }
		//$$ pasteIcon.move(24,312);
		try {
			aboutIcon.setImageLocation(new URL(simulate.Default.IMAGE_DIRECTORY+"about.gif"));
		}
		catch (java.net.MalformedURLException error) { }
        
        //this.setFloatable(false);
        this.setBounds(0, 0, 200, 56);
        
        /**
         * Add buttons, tooltips, actionListeners, and mnemonics to the JToolBar
         */
        newButton.setDefaultCapable(false);
		newButton.setToolTipText("Create a new document");
		newButton.setMnemonic((int)'N');
		newButton.addActionListener(FileActions.NEW_SIMULATION);
		this.add(newButton);
		newButton.setBounds(16,4,51,27);
		openButton.setDefaultCapable(false);
		openButton.setToolTipText("Open an existing document");
		openButton.setMnemonic((int)'O');
		openButton.addActionListener(FileActions.OPEN);
		this.add(openButton);
		openButton.setBounds(67,4,51,27);
		saveButton.setDefaultCapable(false);
		saveButton.setToolTipText("Save the active document");
		saveButton.setMnemonic((int)'S');
		saveButton.addActionListener(FileActions.SEREDIT);
		this.add(saveButton);
		saveButton.setBounds(118,4,51,27);
		this.add(JToolBarSeparator1);
		JToolBarSeparator1.setBounds(169,2,10,5);
		cutButton.setDefaultCapable(false);
		cutButton.setToolTipText("Cut the selection and put it on the Clipboard");
		cutButton.setMnemonic((int)'T');
		cutButton.addActionListener(ButtonActions.CUT);
		this.add(cutButton);
		cutButton.setBounds(179,4,51,27);
		copyButton.setDefaultCapable(false);
		copyButton.setToolTipText("Copy the selection and put it on the Clipboard");
		copyButton.setMnemonic((int)'C');
		copyButton.addActionListener(ButtonActions.COPY);
		this.add(copyButton);
		copyButton.setBounds(230,4,51,27);
		pasteButton.setDefaultCapable(false);
		pasteButton.setToolTipText("Insert Clipboard contents");
		pasteButton.setMnemonic((int)'P');
		pasteButton.addActionListener(ButtonActions.PASTE);
		this.add(pasteButton);
		pasteButton.setBounds(281,4,51,27);
		this.add(JToolBarSeparator2);
		JToolBarSeparator2.setBounds(332,2,10,5);
		aboutButton.setDefaultCapable(false);
		aboutButton.setToolTipText("Display program information, version number and copyright");
		aboutButton.setMnemonic((int)'A');
		aboutButton.addActionListener(HelpActions.ABOUT);
		this.add(aboutButton);
		aboutButton.setBounds(342,4,51,27);

        /**
         * Add icons to each button instance
         */
        saveButton.setIcon(saveIcon);
		newButton.setIcon(newIcon);
		openButton.setIcon(openIcon);
		aboutButton.setIcon(aboutIcon);
		pasteButton.setIcon(pasteIcon);
		cutButton.setIcon(cutIcon);
		copyButton.setIcon(copyIcon);
	}// end of EtomicaToolBar constructor
}// end of EtomicaToolBar class