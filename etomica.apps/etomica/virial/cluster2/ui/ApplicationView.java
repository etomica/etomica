package etomica.virial.cluster2.ui;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;
import javax.swing.WindowConstants;

import com.jgoodies.looks.LookUtils;

public class ApplicationView extends JFrame {

  private static final long     serialVersionUID = -4087118446258055278L;
  // Main Window Preferred Dimensions
  public static final Dimension MWPD_LARGER      = ApplicationUI.MWP_DIM3;
  public static final Dimension MWPD_SMALLER     = ApplicationUI.MWP_DIM2;
  // Main Window Maximized?
  public static final Boolean   MW_MAXIMIZED     = Boolean.TRUE;
  // Main Window Icon
  public static final String    MW_ICON          = "home.gif";
  // Main Window Title string
  public static final String    MW_TITLE         = "Cluster Diagram Manager Application";
  // Main Window About string choice
  public static final String    MW_ABOUT_TEXT    = "UI Prototype for the Cluster Diagram Manager Application\n\n";
  // Message to indicate a feature has not been implemented yet
  public static final String    NOT_IMPLEMENTED  = "This feature has not been implemented yet.";

  private ApplicationView() {

    createContents();
    configureFrame();
  }

  /**************************************************************
   * Content Creation:
   **************************************************************/
  private void createContents() {

    setContentPane(buildContentPane());
    setJMenuBar(new MenuBarView().buildMenuBar(this));
  }

  /**
   * Creates the client area of the main window. A toolbar is placed at the top
   * of the client area and a main panel occupies the rest of the client area.
   */
  private JComponent buildContentPane() {

    JPanel panel = new JPanel(new BorderLayout());
    panel.add(ToolBarView.mainToolbar(this), BorderLayout.NORTH);
    panel.add(MainPanelView.build(), BorderLayout.CENTER);
    return panel;
  }

  /**************************************************************
   * Frame Specific Configuration
   **************************************************************/
  private void configureFrame() {

    setTitle(MW_TITLE);
    setIconImage(ApplicationUI.readImageIcon(MW_ICON).getImage());
    setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
    setVisible(true);
    setDefaultSizeAndPosition();
  }

  private Dimension getPreferredDefaultSize() {

    return LookUtils.IS_LOW_RESOLUTION ? MWPD_SMALLER : MWPD_LARGER;
  }

  private void setDefaultSizeAndPosition() {

    // set the frame size as the preferred default size; if the window starts
    // maximized, it will later restore to the preferred default size;
    setSize(getPreferredDefaultSize());
    // center the frame on the screen; if the window starts maximized, it will
    // later restore to this position;
    Dimension paneSize = getSize();
    Dimension screenSize = getToolkit().getScreenSize();
    setLocation((screenSize.width - paneSize.width) / 2,
        (screenSize.height - paneSize.height) / 2);
    // maximize the frame if necessary
    if (MW_MAXIMIZED.booleanValue()) {
      setExtendedState(JFrame.MAXIMIZED_BOTH);
    }
  }

  /**************************************************************
   * Public Listener Inner Classes
   **************************************************************/
  /**
   * Creates an ActionListener that opens opens up the help.
   */
  final class HelpActionListener implements ActionListener {

    public void actionPerformed(ActionEvent e) {

      JOptionPane.showMessageDialog(ApplicationView.this,
          ApplicationView.NOT_IMPLEMENTED);
    }
  }

  /**
   * Creates an ActionListener that opens the about dialog.
   */
  final class AboutActionListener implements ActionListener {

    public void actionPerformed(ActionEvent e) {

      JOptionPane.showMessageDialog(ApplicationView.this,
          ApplicationView.MW_ABOUT_TEXT);
    }
  }

  /**
   * Creates an ActionListener that opens the file chooser dialog.
   */
  final class OpenFileActionListener implements ActionListener {

    public void actionPerformed(ActionEvent e) {

      new JFileChooser().showOpenDialog(ApplicationView.this);
    }
  }

  /**************************************************************
   * Static Entry Point to Run the Application
   **************************************************************/
  public static void showMainWindow() {

    ApplicationUI.configure();
    SwingUtilities.invokeLater(new Runnable() {

      public void run() {

        new ApplicationView();
      }
    });
  }
}