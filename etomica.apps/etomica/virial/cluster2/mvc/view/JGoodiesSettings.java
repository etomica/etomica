package etomica.virial.cluster2.mvc.view;

import javax.swing.LookAndFeel;

import com.jgoodies.looks.BorderStyle;
import com.jgoodies.looks.HeaderStyle;
import com.jgoodies.looks.plastic.PlasticLookAndFeel;
import com.jgoodies.looks.plastic.PlasticTheme;
import com.jgoodies.looks.plastic.PlasticXPLookAndFeel;

/**
 * Describes most of the optional settings of the JGoodies Looks. Used by the
 * <code>DemoFrame</code> to configure the UI.
 * 
 * @author Karsten Lentzsch
 * @version $Revision$
 * @see com.jgoodies.looks.BorderStyle
 * @see com.jgoodies.looks.HeaderStyle
 * @see com.jgoodies.looks.Options
 */
public final class JGoodiesSettings {

  private LookAndFeel  selectedLookAndFeel;
  private PlasticTheme selectedTheme;
  private boolean      useNarrowButtons;
  private boolean      tabIconsEnabled;
  private String       plasticTabStyle;
  private boolean      plasticHighContrastFocusEnabled;
  private Boolean      popupDropShadowEnabled;
  private HeaderStyle  menuBarHeaderStyle;
  private BorderStyle  menuBarPlasticBorderStyle;
  private BorderStyle  menuBarWindowsBorderStyle;
  private Boolean      menuBar3DHint;
  private HeaderStyle  toolBarHeaderStyle;
  private BorderStyle  toolBarPlasticBorderStyle;
  private BorderStyle  toolBarWindowsBorderStyle;
  private Boolean      toolBar3DHint;

  private JGoodiesSettings() {

    // Override default constructor; prevents instantiation.
  }

  public static JGoodiesSettings createDefault() {

    JGoodiesSettings settings = new JGoodiesSettings();
    settings.setSelectedLookAndFeel(new PlasticXPLookAndFeel());
    settings.setSelectedTheme(PlasticLookAndFeel.createMyDefaultTheme());
    settings.setUseNarrowButtons(true);
    settings.setTabIconsEnabled(true);
    settings.setPlasticTabStyle(PlasticLookAndFeel.TAB_STYLE_DEFAULT_VALUE);
    settings.setPlasticHighContrastFocusEnabled(false);
    settings.setPopupDropShadowEnabled(null);
    settings.setMenuBarHeaderStyle(null);
    settings.setMenuBarPlasticBorderStyle(null);
    settings.setMenuBarWindowsBorderStyle(null);
    settings.setMenuBar3DHint(null);
    settings.setToolBarHeaderStyle(null);
    settings.setToolBarPlasticBorderStyle(null);
    settings.setToolBarWindowsBorderStyle(null);
    settings.setToolBar3DHint(null);
    return settings;
  }

  public Boolean getMenuBar3DHint() {

    return menuBar3DHint;
  }

  public void setMenuBar3DHint(Boolean menuBar3DHint) {

    this.menuBar3DHint = menuBar3DHint;
  }

  public HeaderStyle getMenuBarHeaderStyle() {

    return menuBarHeaderStyle;
  }

  public void setMenuBarHeaderStyle(HeaderStyle menuBarHeaderStyle) {

    this.menuBarHeaderStyle = menuBarHeaderStyle;
  }

  public BorderStyle getMenuBarPlasticBorderStyle() {

    return menuBarPlasticBorderStyle;
  }

  public void setMenuBarPlasticBorderStyle(BorderStyle menuBarPlasticBorderStyle) {

    this.menuBarPlasticBorderStyle = menuBarPlasticBorderStyle;
  }

  public BorderStyle getMenuBarWindowsBorderStyle() {

    return menuBarWindowsBorderStyle;
  }

  public void setMenuBarWindowsBorderStyle(BorderStyle menuBarWindowsBorderStyle) {

    this.menuBarWindowsBorderStyle = menuBarWindowsBorderStyle;
  }

  public Boolean isPopupDropShadowEnabled() {

    return popupDropShadowEnabled;
  }

  public void setPopupDropShadowEnabled(Boolean popupDropShadowEnabled) {

    this.popupDropShadowEnabled = popupDropShadowEnabled;
  }

  public boolean isPlasticHighContrastFocusEnabled() {

    return plasticHighContrastFocusEnabled;
  }

  public void setPlasticHighContrastFocusEnabled(
      boolean plasticHighContrastFocusEnabled) {

    this.plasticHighContrastFocusEnabled = plasticHighContrastFocusEnabled;
  }

  public String getPlasticTabStyle() {

    return plasticTabStyle;
  }

  public void setPlasticTabStyle(String plasticTabStyle) {

    this.plasticTabStyle = plasticTabStyle;
  }

  public LookAndFeel getSelectedLookAndFeel() {

    return selectedLookAndFeel;
  }

  public void setSelectedLookAndFeel(LookAndFeel selectedLookAndFeel) {

    this.selectedLookAndFeel = selectedLookAndFeel;
  }

  public void setSelectedLookAndFeel(String selectedLookAndFeelClassName) {

    try {
      Class theClass = Class.forName(selectedLookAndFeelClassName);
      setSelectedLookAndFeel((LookAndFeel) theClass.newInstance());
    }
    catch (Exception e) {
      System.out.println("Can't instantiate " + selectedLookAndFeelClassName);
      e.printStackTrace();
    }
  }

  public PlasticTheme getSelectedTheme() {

    return selectedTheme;
  }

  public void setSelectedTheme(PlasticTheme selectedTheme) {

    this.selectedTheme = selectedTheme;
  }

  public boolean isTabIconsEnabled() {

    return tabIconsEnabled;
  }

  public void setTabIconsEnabled(boolean tabIconsEnabled) {

    this.tabIconsEnabled = tabIconsEnabled;
  }

  public Boolean getToolBar3DHint() {

    return toolBar3DHint;
  }

  public void setToolBar3DHint(Boolean toolBar3DHint) {

    this.toolBar3DHint = toolBar3DHint;
  }

  public HeaderStyle getToolBarHeaderStyle() {

    return toolBarHeaderStyle;
  }

  public void setToolBarHeaderStyle(HeaderStyle toolBarHeaderStyle) {

    this.toolBarHeaderStyle = toolBarHeaderStyle;
  }

  public BorderStyle getToolBarPlasticBorderStyle() {

    return toolBarPlasticBorderStyle;
  }

  public void setToolBarPlasticBorderStyle(BorderStyle toolBarPlasticBorderStyle) {

    this.toolBarPlasticBorderStyle = toolBarPlasticBorderStyle;
  }

  public BorderStyle getToolBarWindowsBorderStyle() {

    return toolBarWindowsBorderStyle;
  }

  public void setToolBarWindowsBorderStyle(BorderStyle toolBarWindowsBorderStyle) {

    this.toolBarWindowsBorderStyle = toolBarWindowsBorderStyle;
  }

  public boolean isUseNarrowButtons() {

    return useNarrowButtons;
  }

  public void setUseNarrowButtons(boolean useNarrowButtons) {

    this.useNarrowButtons = useNarrowButtons;
  }
}