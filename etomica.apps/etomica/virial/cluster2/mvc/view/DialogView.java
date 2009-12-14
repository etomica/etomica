package etomica.virial.cluster2.mvc.view;

import java.awt.Dimension;
import java.awt.Frame;

import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.WindowConstants;

import com.jgoodies.looks.LookUtils;

import etomica.virial.cluster2.mvc.State;
import etomica.virial.cluster2.mvc.View;
import etomica.virial.cluster2.mvc.ViewResponse;

public abstract class DialogView extends JDialog implements View {

  private static final long serialVersionUID = 1988692491467886360L;
  // Dialog's Preferred Dimensions
// public static final Dimension DP_DIM1 = new Dimension(560, 420);
// public static final Dimension DP_DIM2 = new Dimension(640, 480);
// public static final Dimension DP_DIM3 = new Dimension(720, 540);
// public static final Dimension DP_DIM4 = new Dimension(800, 600);
  public static final Dimension DP_DIM1 = new Dimension(720, 540);
  public static final Dimension DP_DIM2 = new Dimension(800, 600);
  public static final Dimension DP_DIM3 = new Dimension(880, 660);
  public static final Dimension DP_DIM4 = new Dimension(960, 720);

  public DialogView(final Frame owner, final String dialogTitle) {

    super(owner, dialogTitle);
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

    return DP_DIM3;
  }

  // Determines this dialogs's preferred size in low resolution
  protected Dimension getPreferredLowResSize() {

    return DP_DIM2;
  }

  public void display() {

    setVisible(true);
  }

  public void configure(State state) {

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