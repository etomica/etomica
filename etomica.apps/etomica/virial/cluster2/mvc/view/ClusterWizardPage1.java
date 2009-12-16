package etomica.virial.cluster2.mvc.view;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import com.jgoodies.forms.builder.DefaultFormBuilder;
import com.jgoodies.forms.layout.CellConstraints;
import com.jgoodies.forms.layout.FormLayout;

public class ClusterWizardPage1 extends ClusterWizardPageTemplate {

  private JComponent edClusterName;
  private JSpinner spTotalNodes;
  private JSpinner spRootNodes;
  private JSpinner spFieldNodes;
  private JCheckBox ckIsomorphFree;
  private JComboBox cbColor;
  private JComboBox cbMonochromatic;
  private boolean spinnerChanging = false;

  @Override
  public void attachDone() {

    ((JButton) getData().getProperty(ClusterWizard.KEY_HELP_BUTTON)).setEnabled(false);
    ((JButton) getData().getProperty(ClusterWizard.KEY_BACK_BUTTON)).setEnabled(false);
    ((JButton) getData().getProperty(ClusterWizard.KEY_FINISH_BUTTON)).setEnabled(false);
  }

  @Override
  public void detachDone() {

    ((JButton) getData().getProperty(ClusterWizard.KEY_HELP_BUTTON)).setEnabled(true);
    ((JButton) getData().getProperty(ClusterWizard.KEY_BACK_BUTTON)).setEnabled(true);
    ((JButton) getData().getProperty(ClusterWizard.KEY_FINISH_BUTTON)).setEnabled(true);
  }

  @Override
  protected int getPageIndex() {

    return 1;
  }

  @Override
  protected JComponent createControls() {

    FormLayout layout = new FormLayout("250dlu", "pref, 10dlu:grow, pref, 20dlu:grow, pref, 10dlu:grow, pref");
    DefaultFormBuilder builder = new DefaultFormBuilder(layout);
    builder.setBorder(new EmptyBorder(0, 0, 0, 0));
    builder.setOpaque(false);
    // section
    builder.addSeparator("Global Properties", new CellConstraints(1, 1));
    builder.add(globalSection(), new CellConstraints(1, 3));
    return builder.getPanel();
  }

  protected JComponent globalSection() {

    FormLayout layout = new FormLayout("right:max(100dlu;pref), 6dlu, 80dlu");
    DefaultFormBuilder builder = new DefaultFormBuilder(layout);
    builder.setBorder(new EmptyBorder(0, 0, 0, 0));
    builder.setOpaque(false);

    edClusterName = createText("NewCluster");
    spTotalNodes = createSpinner(new SpinnerNumberModel(4, 1, 13, 1), true, true);
    spRootNodes = createSpinner(new SpinnerNumberModel(2, 0, 4, 1), true, true);
    spFieldNodes = createSpinner(new SpinnerNumberModel(2, 0, 4, 1), false, false);
    ChangeListener listener = new ChangeListener() {

      public void stateChanged(ChangeEvent e) {

        if (spinnerChanging) {
          return;
        }
        spinnerChanging = true;
        try {
          JSpinner spinner = (JSpinner) e.getSource();
          SpinnerNumberModel totalModel = (SpinnerNumberModel) spTotalNodes.getModel();
          SpinnerNumberModel fieldModel = (SpinnerNumberModel) spFieldNodes.getModel();
          SpinnerNumberModel rootModel = (SpinnerNumberModel) spRootNodes.getModel();
          Integer totalValue = (Integer) totalModel.getValue();
          if (spinner == spTotalNodes) {
            rootModel.setMaximum(totalValue);
            if (totalValue < (Integer) rootModel.getValue()) {
              rootModel.setValue(totalValue);
            }
          }
          fieldModel.setValue(totalValue - (Integer) rootModel.getValue());
          fieldModel.setMaximum(totalModel.getMaximum());
        }
        finally {
          spinnerChanging = false;
        }
      }

    };
    spTotalNodes.addChangeListener(listener);
    spRootNodes.addChangeListener(listener);
    ckIsomorphFree = createCheckBox("isomorph-free cluster");
    cbMonochromatic = createComboBox(new String[] { "monochromatic", "multi-colored" }, true);
    cbColor = createComboBox(new String[] { "black", "blue", "red" }, true);
    builder.append("Cluster &Name:", buildGrid(edClusterName));
    builder.append("&Total Nodes:", buildGrid(spTotalNodes));
    builder.append("&Root Nodes:", buildGrid(spRootNodes));
    builder.append("Field Nodes:", buildGrid(spFieldNodes));
    builder.append("Cluster &Colors:", buildGrid(cbMonochromatic));
    builder.append("&Default Node Color:", buildGrid(cbColor));
    builder.append("&Isomorphism:", buildButtonRow(ckIsomorphFree, true, true, false, false));
    return builder.getPanel();
  }

  @Override
  protected String getTitle() {

    String title = "Provide a name for the cluster as well as the node, color and isomorphism specifications. ";
    title += "If the cluster is monochromatic, field and root nodes will all have the same color (the default color). ";
    title += "If the cluster is multi-colored, you will be asked to define the colors of the root and field nodes before the cluster is generated.";
    return title;
  }
}