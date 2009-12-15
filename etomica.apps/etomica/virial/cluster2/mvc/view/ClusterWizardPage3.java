package etomica.virial.cluster2.mvc.view;

import java.awt.Component;
import javax.swing.BorderFactory;
import javax.swing.DefaultCellEditor;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.SwingConstants;
import javax.swing.border.EmptyBorder;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.TableColumn;

import com.jgoodies.forms.builder.DefaultFormBuilder;
import com.jgoodies.forms.layout.CellConstraints;
import com.jgoodies.forms.layout.FormLayout;
import com.jgoodies.looks.LookUtils;
import com.jgoodies.uif_lite.component.Factory;

public class ClusterWizardPage3 extends ClusterWizardPageTemplate {

  @Override
  public void attachDone() {

    ((JButton) getData().getProperty(ClusterWizard.KEY_HELP_BUTTON)).setEnabled(false);
    ((JButton) getData().getProperty(ClusterWizard.KEY_FINISH_BUTTON)).setEnabled(false);
  }

  @Override
  public void detachDone() {

    ((JButton) getData().getProperty(ClusterWizard.KEY_HELP_BUTTON)).setEnabled(true);
    ((JButton) getData().getProperty(ClusterWizard.KEY_FINISH_BUTTON)).setEnabled(true);
  }

  @Override
  protected int getPageIndex() {

    return 3;
  }

  @Override
  protected JComponent createControls() {

    FormLayout layout = new FormLayout("250dlu", "pref, 10dlu:grow, pref");
    DefaultFormBuilder builder = new DefaultFormBuilder(layout);
    builder.setBorder(new EmptyBorder(0, 0, 0, 0));
    builder.setOpaque(false);
    // section
    builder.addSeparator("Color Specifications", new CellConstraints(1, 1));
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

    String[] columnNames = { "id", "name", "color" };
    Object[][] data = { { "1", "root node #1", "Black" }, { "2", "root node #2", "Black" },
        { "3", "field node #1", "Yellow" }, { "4", "field node #2", "Red" }, { "5", "field node #3", "Red" },
        { "6", "field node #4", "Blue" } };

    JTable table = new JTable(data, columnNames);
    // table.setBorder(new EmptyBorder(2, 2, 2, 2));
    table.setOpaque(false);

    table.getColumnModel().getColumn(0).setPreferredWidth(30);
    table.getColumnModel().getColumn(1).setPreferredWidth(150);
    table.getColumnModel().getColumn(2).setPreferredWidth(150);
    int tableFontSize = table.getFont().getSize();
    int minimumRowHeight = tableFontSize + 10;
    int defaultRowHeight = LookUtils.IS_LOW_RESOLUTION ? 24 : 32;
    table.setRowHeight(Math.max(minimumRowHeight, defaultRowHeight));

    TableColumn editableColumn = table.getColumnModel().getColumn(2);
    JComboBox comboBox = createComboBox(new String[] { "Black", "Blue", "Red", "Yellow" }, true);
    editableColumn.setCellEditor(new DefaultCellEditor(comboBox));

    // let's add an internal left margin on text cells
    table.setDefaultRenderer(Object.class, new DefaultTableCellRenderer() {

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

    // let's add an internal right margin on numeric cells
    table.setDefaultRenderer(Integer.class, new DefaultTableCellRenderer() {

      private static final long serialVersionUID = -614070586988865358L;

      @Override
      public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected,
          boolean hasFocus, int row, int column) {

        super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, column);
        setBorder(BorderFactory
            .createCompoundBorder(getBorder(), BorderFactory.createEmptyBorder(0, 0, 0, 8)));
        return this;
      }
    });

    // associate the numeric renderer with the first column, and make sure we align right
    TableColumn numericColumn = table.getColumnModel().getColumn(0);
    numericColumn.setCellRenderer(table.getDefaultRenderer(Integer.class));
    ((JLabel) numericColumn.getCellRenderer()).setHorizontalAlignment(SwingConstants.RIGHT);

    JScrollPane pane = Factory.createStrippedScrollPane(table);
    pane.setOpaque(false);
    pane.getViewport().setOpaque(false);

    return pane;
  }

  @Override
  protected String getTitle() {

    String title = "Specify the node color properties for the cluster. Note that if you specified a monochromatic cluster ";
    title += "you should change the node colors by using the default color combo box on the first page of this wizard.";
    return title;
  }
}