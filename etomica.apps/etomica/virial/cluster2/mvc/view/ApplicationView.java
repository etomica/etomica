/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.mvc.view;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JMenuBar;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

public class ApplicationView extends FrameView {

  private static final long serialVersionUID = -4087118446258055278L;
  // Main Window Icon
  public static final String MW_ICON = "home.gif";
  // Main Window Icon
  public static final boolean MW_MAXIMIZED = true;
  // Main Window Title string
  public static final String MW_TITLE = "Cluster Diagram Manager Application";
  // Main Window About string choice
  public static final String MW_ABOUT_TEXT = "UI Prototype for the Cluster Diagram Manager Application\n\n";
  // Message to indicate a feature has not been implemented yet
  public static final String NOT_IMPLEMENTED = "This feature has not been implemented yet.";

  public ApplicationView() {

    super(MW_TITLE);
  }

  /**
   * Creates the client area of the main window. A toolbar is placed at the top
   * of the client area and a main panel occupies the rest of the client area.
   */
  @Override
  protected JComponent buildContentPane() {

    JPanel panel = new JPanel(new BorderLayout());
    panel.add(ToolBarView.mainToolbar(this), BorderLayout.NORTH);
    panel.add(MainPanelView.build(), BorderLayout.CENTER);
    return panel;
  }

  // Creates the main menu of the frame.
  @Override
  protected JMenuBar buildMainMenu() {

    return (new MenuBarView()).buildMenuBar(this);
  }

  // Returns the resource name of this frame's icon
  @Override
  protected String getIconName() {

    return MW_ICON;
  }

  // Determines this frame's preferred size in high resolution
  @Override
  protected Dimension getPreferredHighResSize() {

    return FP_DIM4;
  }

  // Determines this frame's preferred size in low resolution
  @Override
  protected Dimension getPreferredLowResSize() {

    return FP_DIM3;
  }

  // Determines if this frame should display maximized
  @Override
  protected boolean isMaximized() {

    return MW_MAXIMIZED;
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
}