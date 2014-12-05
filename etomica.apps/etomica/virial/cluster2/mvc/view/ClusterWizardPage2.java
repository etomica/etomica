/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.mvc.view;

import java.awt.Component;

import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.*;

import com.jgoodies.forms.builder.DefaultFormBuilder;
import com.jgoodies.forms.layout.CellConstraints;
import com.jgoodies.forms.layout.FormLayout;

import etomica.virial.cluster2.mvc.WizardController;

import static etomica.virial.cluster2.mvc.view.ClusterWizardState.*;

public class ClusterWizardPage2 extends ClusterWizardPageTemplate {

  public ClusterWizardPage2(WizardController controller) {

    super(controller);
  }

  private JRadioButton rbConnectivityReeHoover;
  private JRadioButton rbConnectivityBiconnected;
  private JRadioButton rbConnectivityConnected;
  private JRadioButton rbConnectivityAny;
  private JCheckBox ckNodalPoint;
  private JCheckBox ckArticulationPoint;
  private JCheckBox ckArticulationPair;

  @Override
  public void attachDone() {

    ((JButton) getController().getState().getProperty(ClusterWizard.KEY_HELP_BUTTON)).setEnabled(false);
    ((JButton) getController().getState().getProperty(ClusterWizard.KEY_FINISH_BUTTON)).setEnabled(false);
    super.attachDone();
  }

  @Override
  public void detachDone() {

    ((JButton) getController().getState().getProperty(ClusterWizard.KEY_HELP_BUTTON)).setEnabled(true);
    ((JButton) getController().getState().getProperty(ClusterWizard.KEY_FINISH_BUTTON)).setEnabled(true);
    super.detachDone();
  }

  @Override
  public int getPageId() {

    return 2;
  }

  @Override
  protected JComponent createControls() {

    FormLayout layout = new FormLayout("250dlu", "pref, 10dlu:grow, pref, 20dlu:grow, pref, 10dlu:grow, pref");
    DefaultFormBuilder builder = new DefaultFormBuilder(layout);
    builder.setBorder(new EmptyBorder(0, 0, 0, 0));
    builder.setOpaque(false);
    // section
    builder.addSeparator("Connectivity Class", new CellConstraints(1, 1));
    builder.add(connectivitySection(), new CellConstraints(1, 3));
    // section
    builder.addSeparator("Connectivity Filters", new CellConstraints(1, 5));
    builder.add(filtersSection(), new CellConstraints(1, 7));
    return builder.getPanel();
  }

  protected Component connectivitySection() {

    FormLayout layout = new FormLayout("30dlu, 190dlu", "pref, pref, pref, pref");
    DefaultFormBuilder builder = new DefaultFormBuilder(layout);
    builder.setBorder(new EmptyBorder(0, 0, 0, 0));
    builder.setOpaque(false);
    rbConnectivityAny = createRadioButton("all clusters");
    rbConnectivityConnected = createRadioButton("connected clusters only");
    rbConnectivityBiconnected = createRadioButton("biconnected clusters only");
    rbConnectivityReeHoover = createRadioButton("Ree-Hoover clusters only");
    ButtonGroup group = new ButtonGroup();
    group.add(rbConnectivityAny);
    group.add(rbConnectivityConnected);
    group.add(rbConnectivityBiconnected);
    group.add(rbConnectivityReeHoover);
    builder.add(buildButtonRow(rbConnectivityAny, true, true, false, false), new CellConstraints(2, 1));
    builder
        .add(buildButtonRow(rbConnectivityConnected, false, true, false, false), new CellConstraints(2, 2));
    builder.add(buildButtonRow(rbConnectivityBiconnected, false, true, false, false), new CellConstraints(2,
        3));
    builder
        .add(buildButtonRow(rbConnectivityReeHoover, false, true, false, false), new CellConstraints(2, 4));
    return builder.getPanel();
  }

  protected Component filtersSection() {

    FormLayout layout = new FormLayout("30dlu, 190dlu", "pref, pref, pref");
    DefaultFormBuilder builder = new DefaultFormBuilder(layout);
    builder.setBorder(new EmptyBorder(0, 0, 0, 0));
    builder.setOpaque(false);
    ckNodalPoint = createCheckBox("exclude clusters with nodal points");
    ckArticulationPoint = createCheckBox("exclude clusters with articulation points");
    ckArticulationPair = createCheckBox("exclude clusters with articulation pairs");
    builder.add(buildButtonRow(ckNodalPoint, false, false, false, false), new CellConstraints(2, 1));
    builder.add(buildButtonRow(ckArticulationPoint, false, false, false, false), new CellConstraints(2, 2));
    builder.add(buildButtonRow(ckArticulationPair, false, false, false, false), new CellConstraints(2, 3));
    ChangeListener listener = new ChangeListener() {

      public void stateChanged(ChangeEvent e) {

        ckNodalPoint.setEnabled(rbConnectivityConnected.isSelected());
        ckArticulationPoint.setEnabled(rbConnectivityConnected.isSelected());
        ckArticulationPair.setEnabled(rbConnectivityConnected.isSelected()
            || rbConnectivityBiconnected.isSelected());
      }
    };
    rbConnectivityAny.addChangeListener(listener);
    rbConnectivityConnected.addChangeListener(listener);
    rbConnectivityBiconnected.addChangeListener(listener);
    rbConnectivityReeHoover.addChangeListener(listener);
    return builder.getPanel();
  }

  @Override
  protected String getTitle() {

    String title = "Specify the connectivity properties for the cluster. Four main classes of connectivty are supported: ";
    title += "any, connected, biconnected, and Ree-Hoover. Note that connected and biconnected here refer to the graph ";
    title += "theoretic definitions. For connected and biconnected clusters, you may specify two additional filters. ";
    return title;
  }

  @Override
  public void commitChanges() {

    getController().getState().setProperty(KEY_CLASS_ANY, rbConnectivityAny.isSelected());
    getController().getState().setProperty(KEY_CLASS_CONNECTED, rbConnectivityConnected.isSelected());
    getController().getState().setProperty(KEY_CLASS_BICONNECTED, rbConnectivityBiconnected.isSelected());
    getController().getState().setProperty(KEY_CLASS_REEHOOVER, rbConnectivityReeHoover.isSelected());
    getController().getState().setProperty(KEY_EXCLUDE_NODAL_POINTS, ckNodalPoint.isSelected());
    getController().getState().setProperty(KEY_EXCLUDE_ARTICULATION_POINTS, ckArticulationPoint.isSelected());
    getController().getState().setProperty(KEY_EXCLUDE_ARTICULATION_PAIRS, ckArticulationPair.isSelected());
  }
}