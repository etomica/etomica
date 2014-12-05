/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.mvc.view;

import java.awt.Component;
import java.util.ArrayList;
import java.util.List;

import javax.swing.*;
import javax.swing.border.*;
import javax.swing.table.*;

import com.jgoodies.forms.builder.DefaultFormBuilder;
import com.jgoodies.forms.layout.CellConstraints;
import com.jgoodies.forms.layout.FormLayout;
import com.jgoodies.looks.LookUtils;
import com.jgoodies.uif_lite.component.Factory;

import etomica.virial.cluster2.mvc.WizardController;

import static etomica.virial.cluster2.mvc.view.ClusterWizardState.*;

public class ClusterWizardPage4 extends ClusterWizardPageTemplate {

  private JTable colorAssignment;

  public ClusterWizardPage4(WizardController controller) {

    super(controller);
  }

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

    return 4;
  }

  @Override
  protected JComponent createControls() {

    FormLayout layout = new FormLayout("250dlu", "pref, 10dlu:grow, pref");
    DefaultFormBuilder builder = new DefaultFormBuilder(layout);
    builder.setBorder(new EmptyBorder(0, 0, 0, 0));
    builder.setOpaque(false);
    // section
    builder.addSeparator("Color Assignments", new CellConstraints(1, 1));
    builder.add(colorsSection(), new CellConstraints(1, 3));
    return builder.getPanel();
  }

  protected Component colorsSection() {

    FormLayout layout = new FormLayout("30dlu, 190dlu", "160dlu:grow");
    DefaultFormBuilder builder = new DefaultFormBuilder(layout);
    builder.setBorder(new EmptyBorder(0, 0, 0, 0));
    builder.setOpaque(false);
    builder.add(createColorTable(), new CellConstraints(2, 1));
    return builder.getPanel();
  }

  @SuppressWarnings("unchecked")
  protected JComponent createColorTable() {

    String[] columnNames = { "node id", "color" };
    int totalNodes = (Integer) getController().getState().getProperty(KEY_TOTAL_NODES);
    int rootNodes = (Integer) getController().getState().getProperty(KEY_ROOT_NODES);
    List<ColorEntry> mappedColors = (List<ColorEntry>) getController().getState().getProperty(
        KEY_MAPPED_COLORS);
    int colorIndex = 0;
    Object[][] data = new Object[totalNodes][];
    // create the root nodes
    for (int i = 0; i < rootNodes; i++) {
      data[i] = new Object[2];
      data[i][0] = String.format("root node %d", i);
      data[i][1] = mappedColors.get(colorIndex);
      if (colorIndex < (mappedColors.size() - 1)) {
        colorIndex++;
      }
    }
    // create the field nodes
    for (int i = 0; i < totalNodes - rootNodes; i++) {
      data[rootNodes + i] = new Object[2];
      data[rootNodes + i][0] = String.format("field node %d", i);
      data[rootNodes + i][1] = mappedColors.get(colorIndex);
      if (colorIndex < (mappedColors.size() - 1)) {
        colorIndex++;
      }
    }

    colorAssignment = new JTable(data, columnNames);
    // table.setBorder(new EmptyBorder(2, 2, 2, 2));
    colorAssignment.setOpaque(false);

    colorAssignment.getColumnModel().getColumn(0).setPreferredWidth(40);
    colorAssignment.getColumnModel().getColumn(1).setPreferredWidth(160);
    int tableFontSize = colorAssignment.getFont().getSize();
    int minimumRowHeight = tableFontSize + 10;
    int defaultRowHeight = LookUtils.IS_LOW_RESOLUTION ? 24 : 32;
    colorAssignment.setRowHeight(Math.max(minimumRowHeight, defaultRowHeight));

    TableColumn editableColumn = colorAssignment.getColumnModel().getColumn(1);
    editableColumn.setCellRenderer(new ColorTableCellRenderer());
    editableColumn.setCellEditor(new DefaultCellEditor(createColorComboBox(mappedColors)));

    // let's add an internal left margin on text cells
    colorAssignment.setDefaultRenderer(Object.class, new DefaultTableCellRenderer() {

      private static final long serialVersionUID = 2329298432184206803L;

      @Override
      public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected,
          boolean hasFocus, int row, int column) {

        super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, column);
        setBorder(BorderFactory
            .createCompoundBorder(getBorder(), BorderFactory.createEmptyBorder(0, 8, 0, 0)));
        return this;
      }
    });

    JScrollPane pane = Factory.createStrippedScrollPane(colorAssignment);
    pane.setOpaque(false);
    pane.getViewport().setOpaque(false);

    return pane;
  }

  @Override
  protected String getTitle() {

    String title = "Assign the colors of the nodes in the cluster. Note that the colors available for ";
    title += "the assignment are those defined by the mapping on the previous page.";
    return title;
  }

  @Override
  public void loadFromState() {

    // this page already loads its data directly from the state
  }

  @Override
  public void commitChanges() {

    TableModel data = colorAssignment.getModel();
    List<ColorEntry> colors = new ArrayList<ColorEntry>();
    // these are the actual (color id, color) pairs of the mapping
    for (int i = 0; i < data.getRowCount(); i++) {
      ColorEntry entry = (ColorEntry) data.getValueAt(i, 1);
      getController().getState().setProperty(KEY_NODE_COLORS.get(i), entry);
      if (!colors.contains(entry)) {
        colors.add(entry);
      }
    }
    // this is the set of colors in the image of the mapping
    getController().getState().setProperty(KEY_ASSIGNED_COLORS, colors);
  }
}