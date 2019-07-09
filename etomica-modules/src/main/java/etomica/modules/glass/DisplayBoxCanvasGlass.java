/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

public interface DisplayBoxCanvasGlass {
    void setConfigIndex(int idx);

    int getConfigIndex();

    void setDrawDisplacement(boolean doDrawDisplacement);

    boolean getDrawDisplacement();

    void setConfigStorage(ConfigurationStorage configStorage);
    ConfigurationStorage getConfigStorage();

    void setFlipDisplacement(boolean flipDisplacement);

    boolean getFlipDisplacement();

}
