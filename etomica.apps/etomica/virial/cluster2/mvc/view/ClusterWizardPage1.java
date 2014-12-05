/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.mvc.view;

import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.*;

import com.jgoodies.forms.builder.DefaultFormBuilder;
import com.jgoodies.forms.layout.CellConstraints;
import com.jgoodies.forms.layout.FormLayout;

import etomica.virial.cluster2.mvc.WizardController;
import static etomica.virial.cluster2.mvc.view.ClusterWizardState.*;

public class ClusterWizardPage1 extends ClusterWizardPageTemplate {

  private JTextField edClusterName;
  private JSpinner spTotalNodes;
  private JSpinner spRootNodes;
  private JSpinner spFieldNodes;
  private JCheckBox ckIsomorphFree;
  private JComboBox cbColorScheme;
  private boolean spinnerChanging = false;

  public ClusterWizardPage1(WizardController controller) {

    super(controller);
  }

  @Override
  public void attachDone() {

    ((JButton) getController().getState().getProperty(ClusterWizard.KEY_HELP_BUTTON)).setEnabled(false);
    ((JButton) getController().getState().getProperty(ClusterWizard.KEY_BACK_BUTTON)).setEnabled(false);
    ((JButton) getController().getState().getProperty(ClusterWizard.KEY_FINISH_BUTTON)).setEnabled(false);
    super.attachDone();
  }

  @Override
  public void detachDone() {

    ((JButton) getController().getState().getProperty(ClusterWizard.KEY_HELP_BUTTON)).setEnabled(true);
    ((JButton) getController().getState().getProperty(ClusterWizard.KEY_BACK_BUTTON)).setEnabled(true);
    ((JButton) getController().getState().getProperty(ClusterWizard.KEY_FINISH_BUTTON)).setEnabled(true);
    super.detachDone();
  }

  @Override
  public int getPageId() {

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
    cbColorScheme = createComboBox(new String[] { DEFVAL_MONOCHROMATIC, DEFVAL_MULTICOLORED }, true);
    builder.append("Cluster &Name:", buildGrid(edClusterName));
    builder.append("&Total Nodes:", buildGrid(spTotalNodes));
    builder.append("&Root Nodes:", buildGrid(spRootNodes));
    builder.append("Field Nodes:", buildGrid(spFieldNodes));
    builder.append("Colors Sc&heme:", buildGrid(cbColorScheme));
    builder.append("&Isomorphism:", buildButtonRow(ckIsomorphFree, true, true, false, false));
    return builder.getPanel();
  }

  @Override
  protected String getTitle() {

    String title = "Provide a name and some basic specifications for the cluster. Connectivity filters can be ";
    title += "defined on page 2, color mappings and assignments on pages 3 and 4. A summary with the complete ";
    title += "specifications will be provided before the cluster is generated.";
    return title;
  }

  @Override
  public void loadFromState() {

    edClusterName.setText((String) getController().getState().getProperty(KEY_NAME));
    spTotalNodes.setValue(getController().getState().getProperty(KEY_TOTAL_NODES));
    spRootNodes.setValue(getController().getState().getProperty(KEY_ROOT_NODES));
    cbColorScheme.setSelectedItem(getController().getState().getProperty(KEY_COLOR_SCHEME));
    ckIsomorphFree.setSelected((Boolean) getController().getState().getProperty(KEY_ISOMORPH_FREE));
  }

  @Override
  public void commitChanges() {

    getController().getState().setProperty(KEY_NAME, edClusterName.getText());
    getController().getState().setProperty(KEY_TOTAL_NODES, spTotalNodes.getValue());
    getController().getState().setProperty(KEY_ROOT_NODES, spRootNodes.getValue());
    getController().getState().setProperty(KEY_COLOR_SCHEME, cbColorScheme.getSelectedItem());
    getController().getState().setProperty(KEY_ISOMORPH_FREE, ckIsomorphFree.isSelected());
  }
}