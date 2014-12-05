/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.mvc;

public interface WizardPageView extends View {

  public WizardController getController();

  public int getPageId();

  public void attach(String key, Object object);

  public void attachDone();

  public void detach(String key, Object object);

  public void detachDone();

  // sets the object to call back upon valid view completion
  public void setResponseListener(ViewResponseListener listener);
}