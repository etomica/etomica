/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.mvc.view;

import java.awt.*;
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

public class ClusterWizardPage3 extends ClusterWizardPageTemplate {

  JTable colorMapping;

  public ClusterWizardPage3(WizardController controller) {

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

    return 3;
  }

  @Override
  protected JComponent createControls() {

    FormLayout layout = new FormLayout("250dlu", "pref, 10dlu:grow, pref");
    DefaultFormBuilder builder = new DefaultFormBuilder(layout);
    builder.setBorder(new EmptyBorder(0, 0, 0, 0));
    builder.setOpaque(false);
    // section
    builder.addSeparator("Color Mapping", new CellConstraints(1, 1));
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

  protected JComponent createColorTable() {

    String[] columnNames = { "color id", "assignment" };
    int totalColors = 1;
    // if the cluster is multicolored, we should support totalNodes color mappings
    if (!getController().getState().getProperty(KEY_COLOR_SCHEME).equals(DEFVAL_MONOCHROMATIC)) {
      totalColors = (Integer) getController().getState().getProperty(KEY_TOTAL_NODES);
    }
    Object[][] data = new Object[totalColors][];
    for (int i = 0; i < totalColors; i++) {
      data[i] = new Object[2];
      data[i][0] = KEY_COLORS.get(i);
      data[i][1] = getController().getState().getProperty(KEY_COLORS.get(i));
    }

    colorMapping = new JTable(data, columnNames);
    // table.setBorder(new EmptyBorder(2, 2, 2, 2));
    colorMapping.setOpaque(false);

    colorMapping.getColumnModel().getColumn(0).setPreferredWidth(40);
    colorMapping.getColumnModel().getColumn(1).setPreferredWidth(160);
    int tableFontSize = colorMapping.getFont().getSize();
    int minimumRowHeight = tableFontSize + 10;
    int defaultRowHeight = LookUtils.IS_LOW_RESOLUTION ? 24 : 32;
    colorMapping.setRowHeight(Math.max(minimumRowHeight, defaultRowHeight));

    TableColumn editableColumn = colorMapping.getColumnModel().getColumn(1);
    editableColumn.setCellRenderer(new ColorTableCellRenderer());
    editableColumn.setCellEditor(new DefaultCellEditor(createColorComboBox(DEFVAL_COLORS)));

    // let's add an internal left margin on text cells
    colorMapping.setDefaultRenderer(Object.class, new DefaultTableCellRenderer() {

      private static final long serialVersionUID = 7636488079743438884L;

      @Override
      public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected,
          boolean hasFocus, int row, int column) {

        super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, column);
        setBorder(BorderFactory
            .createCompoundBorder(getBorder(), BorderFactory.createEmptyBorder(0, 8, 0, 0)));
        return this;
      }
    });

    JScrollPane pane = Factory.createStrippedScrollPane(colorMapping);
    pane.setOpaque(false);
    pane.getViewport().setOpaque(false);

    return pane;
  }

  @Override
  protected String getTitle() {

    String title = "You may adjust the color map below by selecting all the colors that should be available for ";
    title += "the cluster nodes color assignment on the next page.";
    return title;
  }

  @Override
  public void loadFromState() {

    // this page already loads its data directly from the state
  }

  @Override
  public void commitChanges() {

    TableModel data = colorMapping.getModel();
    List<ColorEntry> colors = new ArrayList<ColorEntry>();
    int colorIndex = 0;
    for (int i = 0; i < (Integer) getController().getState().getProperty(KEY_TOTAL_NODES); i++) {
      ColorEntry entry = (ColorEntry) data.getValueAt(colorIndex, 1);
      getController().getState().setProperty(KEY_COLORS.get(i), entry);
      if (!colors.contains(entry)) {
        colors.add(entry);
      }
      if (getController().getState().getProperty(KEY_COLOR_SCHEME).equals(DEFVAL_MULTICOLORED)) {
        colorIndex++;
      }
    }
    // this is the set of colors in the image of the mapping
    getController().getState().setProperty(KEY_MAPPED_COLORS, colors);
    // this is ALSO the set of default colors in the assignment of the mapping
    getController().getState().setProperty(KEY_ASSIGNED_COLORS, colors);
  }
}