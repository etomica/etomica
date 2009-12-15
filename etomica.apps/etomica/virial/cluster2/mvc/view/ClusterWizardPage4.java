package etomica.virial.cluster2.mvc.view;

import java.awt.Component;
import java.awt.Font;

import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.border.EmptyBorder;

import com.jgoodies.forms.builder.DefaultFormBuilder;
import com.jgoodies.forms.layout.CellConstraints;
import com.jgoodies.forms.layout.FormLayout;

public class ClusterWizardPage4 extends ClusterWizardPageTemplate {

  @Override
  public void attachDone() {

    ((JButton) getData().getProperty(ClusterWizard.KEY_HELP_BUTTON)).setEnabled(false);
    ((JButton) getData().getProperty(ClusterWizard.KEY_NEXT_BUTTON)).setEnabled(false);
  }

  @Override
  public void detachDone() {

    ((JButton) getData().getProperty(ClusterWizard.KEY_HELP_BUTTON)).setEnabled(true);
    ((JButton) getData().getProperty(ClusterWizard.KEY_NEXT_BUTTON)).setEnabled(true);
  }

  @Override
  protected int getPageIndex() {

    return 4;
  }

  @Override
  protected JComponent createControls() {

    FormLayout layout = new FormLayout("250dlu", "pref, 10dlu:grow, pref");
    DefaultFormBuilder builder = new DefaultFormBuilder(layout);
    builder.setBorder(new EmptyBorder(0, 0, 0, 0));
    builder.setOpaque(false);
    // section
    builder.addSeparator("Summary and Generation Plan", new CellConstraints(1, 1));
    builder.add(summarySection(), new CellConstraints(1, 3));
    return builder.getPanel();
  }

  protected Component summarySection() {

    FormLayout layout = new FormLayout("10dlu, 230dlu", "220dlu:grow");
    DefaultFormBuilder builder = new DefaultFormBuilder(layout);
    builder.setBorder(new EmptyBorder(0, 0, 0, 0));
    builder.setOpaque(false);
    builder.add(createSummary(), new CellConstraints(2, 1));
    return builder.getPanel();
  }

  protected JComponent createSummary() {

    JTextArea summary = buildTextArea(getSummaryText(), false, true);
    summary.setOpaque(true);
    summary.setBackground(ApplicationUI.uiSettings.getSelectedTheme().getPrimaryControl());
    summary.setFont(new Font("Monospaced", summary.getFont().getStyle(), summary.getFont().getSize() - 1));
    summary.setCaretPosition(0);
    JScrollPane pane = new JScrollPane(summary);
    pane.setOpaque(false);
    pane.getViewport().setOpaque(false);

    return pane;
  }

  private String getSummaryText() {

    String result = "";
    result += "Cluster Specification Summary\n";
    result += "-----------------------------\n\n";
    result += "1. Global Properties\n\n";
    result += "   Cluster Name...........: @NAME\n";
    result += "   Cluster Name...........: @NAME\n";
    result += "   Total Nodes............: @TNODES\n";
    result += "   Root Nodes.............: @RNODES\n";
    result += "   Field Nodes............: @FNODES\n";
    result += "   Cluster Colors.........: @COLORS\n";
    result += "   Default Node Color.....: @DNCOLOR\n";
    result += "   Isomorphism............: @ISOMORPHISM\n\n";
    result += "2. Connectivity Properties\n\n";
    result += "   Connectivity Class.....: @CCLASS\n";
    result += "   Connectivity Filters...: @CFILTERS\n\n";
    result += "3. Color Specifications\n\n";
    result += "   Root Nodes.............: @RCOLORS\n";
    result += "   Field Nodes............: @FCOLORS\n\n\n";
    result += "Cluster Generation Plan\n";
    result += "-----------------------\n\n";
    return result;
  }

  @Override
  protected String getTitle() {

    String title = "Check that the cluster specifications are correct before you start the generation. ";
    title += "Please go over the generation plan to have a rough idea of how the generation is setup ";
    title += "and, therefore, how long it may take to complete.";
    return title;
  }
}