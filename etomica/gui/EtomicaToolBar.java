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

package etomica.gui;

import javax.swing.JButton;
import java.net.URL;
import javax.swing.JToolBar.Separator;
import javax.swing.ImageIcon;

public class EtomicaToolBar extends javax.swing.JToolBar {
    
    /**
     * Icons displayed on editing buttons to distinguish them
     */
    static final ImageIcon newIcon = new ImageIcon();
	static final ImageIcon openIcon = new ImageIcon();
	static final ImageIcon saveIcon = new ImageIcon();
	static final ImageIcon cutIcon = new ImageIcon();
	static final ImageIcon copyIcon = new ImageIcon();
	static final ImageIcon pasteIcon = new ImageIcon();
	static final ImageIcon aboutIcon = new ImageIcon();
	static final ImageIcon startIcon = new ImageIcon();
	static final ImageIcon stopIcon = new ImageIcon();
	static final ImageIcon pauseIcon = new ImageIcon();
	

    /**
     * Editing buttons used on the toolbar along with separators to break them up
     */
	public static final JButton newButton = new JButton();
	public static final JButton openButton = new JButton();
	public static final JButton saveButton = new JButton();
	public static final Separator separator1 = new Separator();
	public static final JButton cutButton = new JButton();
	public static final JButton copyButton = new JButton();
	public static final JButton pasteButton = new JButton();
	public static final Separator separator2 = new Separator();
	public static final JButton aboutButton = new JButton();
	public static final Separator separator3 = new Separator();
	public static final JButton startButton = new JButton();
	public static final JButton stopButton = new JButton();
	public static final JButton pauseButton = new JButton();
    
    /**
     * Constructor for the tool bar of the Etomica environment.  It creates all the editing buttons, adds
     * them to the toolbar, gives them icons and tooltips, and finally adds actionlisteners.
     */
    public EtomicaToolBar() {
        
        /**
         * Find and add icons to each button instance
         */
/*        try {
			newIcon.setImageLocation(new URL(etomica.Default.IMAGE_DIRECTORY+"new.gif"));
			openIcon.setImageLocation(new URL(etomica.Default.IMAGE_DIRECTORY+"open.gif"));
			saveIcon.setImageLocation(new URL(etomica.Default.IMAGE_DIRECTORY+"save.gif"));
			cutIcon.setImageLocation(new URL(etomica.Default.IMAGE_DIRECTORY+"cut.gif"));
			copyIcon.setImageLocation(new URL(etomica.Default.IMAGE_DIRECTORY+"copy.gif"));
			pasteIcon.setImageLocation(new URL(etomica.Default.IMAGE_DIRECTORY+"paste.gif"));
			aboutIcon.setImageLocation(new URL(etomica.Default.IMAGE_DIRECTORY+"about.gif"));
			startIcon.setImageLocation(new URL(etomica.Default.IMAGE_DIRECTORY+"start1.gif"));
			stopIcon.setImageLocation(new URL(etomica.Default.IMAGE_DIRECTORY+"stop1.gif"));
			pauseIcon.setImageLocation(new URL(etomica.Default.IMAGE_DIRECTORY+"pause1.gif"));
		}
		catch (java.net.MalformedURLException error) { }
  */      
        //this.setFloatable(false);
        this.setBounds(0, 0, 363, 56);
        
        /**
         * Add buttons, tooltips, actionListeners, and mnemonics to the JToolBar
         */
        newButton.setDefaultCapable(false);
		newButton.setToolTipText("Create a new document");
		newButton.addActionListener(FileActions.NEW_SIMULATION);
		this.add(newButton);
		newButton.setBounds(16,4,51,27);
		openButton.setDefaultCapable(false);
		openButton.setToolTipText("Open an existing document");
		openButton.addActionListener(FileActions.OPEN);
		this.add(openButton);
		openButton.setBounds(67,4,51,27);
		saveButton.setDefaultCapable(false);
		saveButton.setToolTipText("Save the active document");
		saveButton.addActionListener(FileActions.SEREDIT);
		this.add(saveButton);
		saveButton.setBounds(118,4,51,27);
		this.add(separator1);
		separator1.setBounds(169,2,10,5);
		cutButton.setDefaultCapable(false);
		cutButton.setToolTipText("Cut the selection and put it on the Clipboard");
		cutButton.addActionListener(ButtonActions.CUT);
		this.add(cutButton);
		cutButton.setBounds(179,4,51,27);
		copyButton.setDefaultCapable(false);
		copyButton.setToolTipText("Copy the selection and put it on the Clipboard");
		copyButton.addActionListener(ButtonActions.COPY);
		this.add(copyButton);
		copyButton.setBounds(230,4,51,27);
		pasteButton.setDefaultCapable(false);
		pasteButton.setToolTipText("Insert Clipboard contents");
		pasteButton.addActionListener(ButtonActions.PASTE);
		this.add(pasteButton);
		pasteButton.setBounds(281,4,51,27);
		this.add(separator2);
		separator2.setBounds(332,2,10,5);
		aboutButton.setDefaultCapable(false);
		aboutButton.setToolTipText("Display program information, version number and copyright");
		aboutButton.addActionListener(HelpActions.ABOUT);
		this.add(aboutButton);
		aboutButton.setBounds(342,4,51,27);
		this.add(separator3);
		separator3.setBounds(393,2,10,5);
		startButton.setDefaultCapable(false);
		startButton.setToolTipText("Starts the simulation");
		startButton.addActionListener(ControllerActions.START);
		this.add(startButton);
		startButton.setBounds(403,4,51,27);
		stopButton.setDefaultCapable(false);
		stopButton.setToolTipText("Stops the simulation");
		stopButton.addActionListener(ControllerActions.STOP);
		this.add(stopButton);
		stopButton.setBounds(454,4,51,27);
		pauseButton.setDefaultCapable(false);
		pauseButton.setToolTipText("Pauses the simulation");
		pauseButton.addActionListener(ControllerActions.PAUSE);
		this.add(pauseButton);
		pauseButton.setBounds(505,4,51,27);
        pauseButton.setPressedIcon(pauseButton.getDisabledIcon());
        /**
         * Add icons to each button instance
         */
        try {
            saveButton.setIcon(new ImageIcon(new URL(etomica.Default.IMAGE_DIRECTORY+"save.gif")));//saveIcon);
		    newButton.setIcon(new ImageIcon(new URL(etomica.Default.IMAGE_DIRECTORY+"new.gif")));//newIcon);
		    openButton.setIcon(new ImageIcon(new URL(etomica.Default.IMAGE_DIRECTORY+"open.gif")));//openIcon);
		    aboutButton.setIcon(new ImageIcon(new URL(etomica.Default.IMAGE_DIRECTORY+"about.gif")));//aboutIcon);
		    pasteButton.setIcon(new ImageIcon(new URL(etomica.Default.IMAGE_DIRECTORY+"paste.gif")));//pasteIcon);
		    cutButton.setIcon(new ImageIcon(new URL(etomica.Default.IMAGE_DIRECTORY+"cut.gif")));//cutIcon);
		    copyButton.setIcon(new ImageIcon(new URL(etomica.Default.IMAGE_DIRECTORY+"copy.gif")));//copyIcon);
		    startButton.setIcon(new ImageIcon(new URL(etomica.Default.IMAGE_DIRECTORY+"start.gif")));//startIcon);
		    stopButton.setIcon(new ImageIcon(new URL(etomica.Default.IMAGE_DIRECTORY+"stop.gif")));//stopIcon);
		    pauseButton.setIcon(new ImageIcon(new URL(etomica.Default.IMAGE_DIRECTORY+"pause.gif")));//pauseIcon);
		}
		catch (java.net.MalformedURLException me){}
	}// end of EtomicaToolBar constructor
}// end of EtomicaToolBar class