/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.mvc;

public interface WizardState extends State {

  // sets the default state for the page
  public void loadDefaultState(int pageId);

  // determines if the state for the given page has been previously loaded
  public boolean isStateLoaded(int pageId);
}