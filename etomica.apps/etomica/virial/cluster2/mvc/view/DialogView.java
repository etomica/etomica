/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.mvc.view;

import java.awt.Dimension;
import java.awt.Frame;

import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.WindowConstants;

import com.jgoodies.looks.LookUtils;

import etomica.virial.cluster2.mvc.View;

public abstract class DialogView extends JDialog implements View {

  private static final long serialVersionUID = 1988692491467886360L;
  // Dialog's Preferred Dimensions
  public static final Dimension DP_DIM1 = new Dimension(560, 420);
  public static final Dimension DP_DIM2 = new Dimension(640, 480);
  public static final Dimension DP_DIM3 = new Dimension(720, 540);
  public static final Dimension DP_DIM4 = new Dimension(800, 600);
  public static final Dimension DP_DIM5 = new Dimension(880, 660);
  public static final Dimension DP_DIM6 = new Dimension(960, 720);

  public DialogView(final Frame owner, final String dialogTitle) {

    super(owner, dialogTitle);
    // because our DialogView extends JDialog, the background color picked
    // from the theme is not the same as if we had used a JFrame; we use an
    // extra setting to decide whether to override this setting
    if (ApplicationUI.overrideDialogBackground) {
      setBackground(ApplicationUI.uiSettings.getSelectedTheme().getControl());
    }
  }

  public DialogView(final String dialogTitle) {

    this(null, dialogTitle);
  }

  // Creates the client area of the dialog.
  protected abstract JComponent buildContentPane();

  // Returns the resource name of this dialog's icon
  protected String getIconName() {

    return null;
  }

  // Determines this dialogs's preferred size in high resolution
  protected Dimension getPreferredHighResSize() {

    return DP_DIM6;
  }

  // Determines this dialogs's preferred size in low resolution
  protected Dimension getPreferredLowResSize() {

    return DP_DIM4;
  }

  public void display() {

    setVisible(true);
  }

  public void initializeUI() {

    setContentPane(buildContentPane());
    if (getIconName() != null) {
      ((JFrame) getParent()).setIconImage(ApplicationUI.readImageIcon(getIconName()).getImage());
    }
    setSize(getPreferredDefaultSize());
    setPreferredSize(getPreferredDefaultSize());
    Dimension paneSize = getSize();
    Dimension screenSize = getToolkit().getScreenSize();
    setLocation((screenSize.width - paneSize.width) / 2,
        (screenSize.height - paneSize.height) / 2);
    setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
  }

  protected Dimension getPreferredDefaultSize() {

    return LookUtils.IS_LOW_RESOLUTION ? getPreferredLowResSize()
        : getPreferredHighResSize();
  }
}