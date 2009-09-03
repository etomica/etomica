package etomica.virial.cluster2.ui;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.Point;

import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.JSeparator;
import javax.swing.border.EmptyBorder;

import com.jgoodies.forms.builder.ButtonBarBuilder2;
import com.jgoodies.forms.builder.PanelBuilder;
import com.jgoodies.forms.layout.CellConstraints;
import com.jgoodies.forms.layout.FormLayout;

// TODO: parameterize all insets and relative dimensions of width and height of all panels
public class WizardView extends DialogView {

  private static final long serialVersionUID = 4676721929353407232L;
  public static final int WZL_NONE = 0x0000;
  public static final int WZL_LEFT_FIGURE = 0x0001;
  public static final int WZL_TOP_FIGURE = 0x0002;
  public static final String IMG_BACKGROUND = "cluster_4.png";
  private JPanel mainPane = null;
  private JPanel wizardPane = null;
  private JPanel modelPane = null;
  private JPanel figurePane = null;
  private JPanel controlBar = null;
  private JButton helpButton = null;
  private JButton backButton = null;
  private JButton nextButton = null;
  private JButton cancelButton = null;
  private JButton finishButton = null;
  private int wizardLayout = WZL_LEFT_FIGURE;
  private boolean localOpaque = true;
  private int delta = 10;
  private String bkgdImage = "";

  protected WizardView(final Frame owner, final String wizardTitle,
      final String resourceName) {

    super(owner, wizardTitle);
    if ((resourceName == null) && (wizardLayout != WZL_NONE)) {
      bkgdImage = IMG_BACKGROUND;
    }
    else {
      bkgdImage = resourceName;
    }
  }

  public WizardView(final String wizardTitle, final String bkgdImage) {

    this(null, wizardTitle, bkgdImage);
  }

  public WizardView(final String wizardTitle) {

    this(null, wizardTitle, null);
  }

  protected JPanel getWizardPane() {

    return wizardPane;
  }

  protected JPanel getModelPane() {

    return modelPane;
  }

  protected JPanel getMainPane() {

    return mainPane;
  }

  public String getBackgroundImage() {

    return bkgdImage;
  }

  protected JPanel getFigurePane() {

    return figurePane;
  }

  protected JPanel getControlBar() {

    return controlBar;
  }

  protected JButton getHelpButton() {

    return helpButton;
  }

  protected JButton getBackButton() {

    return backButton;
  }

  protected JButton getNextButton() {

    return nextButton;
  }

  protected JButton getCancelButton() {

    return cancelButton;
  }

  protected JButton getFinishButton() {

    return finishButton;
  }

  protected int getWizardLayout() {

    return wizardLayout;
  }

  protected void setWizardLayout(int newLayout) {

    wizardLayout = newLayout;
  }

  @Override
  protected JComponent buildContentPane() {

    buildControlBar();
    buildFigurePane();
    buildModelPane();
    buildWizardPane();
    buildMainPane();
    return mainPane;
  }

  protected void buildControlBar() {

    helpButton = new JButton("Help");
    backButton = new JButton("< Back");
    nextButton = new JButton("Next >");
    cancelButton = new JButton("Cancel");
    finishButton = new JButton("Finish");
    ButtonBarBuilder2 builder = new ButtonBarBuilder2();
    builder.setOpaque(localOpaque);
    builder.setBackground(Color.green);
    builder.setBorder(new EmptyBorder(10, 2, 10, 2));
    builder.addUnrelatedGap();
    builder.addButton(helpButton);
    builder.addUnrelatedGap();
    builder.addGlue();
    builder.addButton(new JButton[] { backButton, nextButton, cancelButton,
        finishButton });
    builder.addUnrelatedGap();
    JPanel buttonBar = builder.getPanel();
// controlBar
    controlBar = new JPanel(new BorderLayout());
    controlBar.setOpaque(localOpaque);
    controlBar.setBackground(Color.gray);
    controlBar.setBorder(new EmptyBorder(2, 10, 2, 10));
    controlBar.add(new JSeparator(JSeparator.HORIZONTAL), BorderLayout.NORTH);
    controlBar.add(buttonBar, BorderLayout.CENTER);
  }

  protected void buildFigurePane() {

    Dimension prefSize;
    if (getWizardLayout() == WZL_TOP_FIGURE) {
      prefSize = new Dimension(getPreferredDefaultSize().width - delta, 120);
    }
    else {
      prefSize = new Dimension(200, 370);
    }
    if (getWizardLayout() != WZL_NONE) {
      figurePane = new ApplicationUI.ImagePanel(getBackgroundImage());
      figurePane.setLayout(new BorderLayout());
      if (getWizardLayout() == WZL_TOP_FIGURE) {
        figurePane.add(new JSeparator(JSeparator.HORIZONTAL),
            BorderLayout.SOUTH);
      }
      else {
        figurePane
            .add(new JSeparator(JSeparator.VERTICAL), BorderLayout.EAST);
      }
    }
    else {
      figurePane = new JPanel();
    }
    figurePane.setOpaque(true);
    figurePane.setBackground(Color.white);
    figurePane.setPreferredSize(prefSize);
  }

  protected void buildModelPane() {

    modelPane = new JPanel();
    modelPane.setOpaque(localOpaque);
    modelPane.setBackground(Color.darkGray);
    if (getWizardLayout() == WZL_TOP_FIGURE) {
      modelPane.setPreferredSize(new Dimension(getPreferredDefaultSize().width
          - delta, 300));
    }
    else {
      modelPane.setPreferredSize(new Dimension(350, 390));
    }
  }

  protected void buildMainPane() {

    FormLayout layout = new FormLayout("1dlu, pref:grow, 1dlu",
        "1dlu, fill:pref:grow, 1dlu, pref, 1dlu");
    PanelBuilder builder = new PanelBuilder(layout);
    builder.setOpaque(localOpaque);
    builder.setBackground(Color.red);
    CellConstraints cc = new CellConstraints();
    builder.add(wizardPane, cc.xy(2, 2));
    builder.add(controlBar, cc.xy(2, 4));
    mainPane = builder.getPanel();
  }

  protected void buildWizardPane() {

    String layoutCols = null;
    String layoutRows = null;
    Point ptFigure = null;
    Point ptModel = null;
    switch (getWizardLayout()) {
      case WZL_LEFT_FIGURE: {
        layoutCols = "pref, 1dlu, pref:grow";
        layoutRows = "fill:pref:grow";
        ptFigure = new Point(1, 1);
        ptModel = new Point(3, 1);
        break;
      }
      case WZL_TOP_FIGURE: {
        layoutCols = "fill:pref:grow";
        layoutRows = "pref, 1dlu, fill:pref:grow";
        ptFigure = new Point(1, 1);
        ptModel = new Point(1, 3);
        break;
      }
      case WZL_NONE: {
        layoutCols = "fill:pref:grow";
        layoutRows = "fill:pref:grow";
        ptFigure = null;
        ptModel = new Point(1, 1);
        break;
      }
    }
    FormLayout layout = new FormLayout(layoutCols, layoutRows);
    PanelBuilder builder = new PanelBuilder(layout);
    builder.setOpaque(localOpaque);
    builder.setBackground(Color.yellow);
    CellConstraints cc = new CellConstraints();
    if ((getFigurePane() != null) && (ptFigure != null)) {
      builder.add(figurePane, cc.xy(ptFigure.x, ptFigure.y));
    }
    builder.add(modelPane, cc.xy(ptModel.x, ptModel.y));
    wizardPane = builder.getPanel();
  }

  @Override
  public void activate() {

    super.activate();
// setResizable(false);
  }

  /**************************************************************
   * Demo the Wizard
   **************************************************************/
  public static void main(String[] args) {

    ViewFactory.showView(ViewFactory.VN_WIZARD);
  }
}