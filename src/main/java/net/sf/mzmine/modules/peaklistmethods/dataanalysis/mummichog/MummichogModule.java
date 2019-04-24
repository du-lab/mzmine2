/*
 * Copyright (C) 2018 Du-Lab Team <dulab.binf@gmail.com>
 *
 * This program is free software; you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
 * even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if
 * not, write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.
 */

package net.sf.mzmine.modules.peaklistmethods.dataanalysis.mummichog;

import java.io.File;
import java.util.*;
import javax.annotation.Nonnull;
import net.sf.mzmine.datamodel.*;
import net.sf.mzmine.modules.MZmineModuleCategory;
import net.sf.mzmine.modules.MZmineRunnableModule;
import net.sf.mzmine.parameters.ParameterSet;
import net.sf.mzmine.taskcontrol.Task;
import net.sf.mzmine.util.ExitCode;

public class MummichogModule implements MZmineRunnableModule {

  private static final String MODULE_NAME = "Mummichog calculation";
  private static final String MODULE_DESCRIPTION =
      "Builds an Activity Network after doing pathway and modular analysis";

  @Override
  public @Nonnull String getName() {
    return MODULE_NAME;
  }

  @Override
  public @Nonnull String getDescription() {
    return MODULE_DESCRIPTION;
  }

  @Override
  @Nonnull
  public ExitCode runModule(@Nonnull MZmineProject project, @Nonnull ParameterSet parameters,
      @Nonnull Collection<Task> tasks) {

    double cutoff = parameters.getParameter(MummichogParameters.cutoff).getValue();
    String network = parameters.getParameter(MummichogParameters.network).getValue();
    boolean force_primary_ion =
        parameters.getParameter(MummichogParameters.force_primary_ion).getValue();
    String modeling = parameters.getParameter(MummichogParameters.modeling).getValue();
    File output = parameters.getParameter(MummichogParameters.output).getValue();

    PeakList[] peakLists =
        parameters.getParameter(MummichogParameters.peakLists).getValue().getMatchingPeakLists();

    for (PeakList peakList : peakLists) {

      tasks.add(new MummichogTask(peakList.getRows(), cutoff, network, force_primary_ion, modeling,
          output));
    }

    return ExitCode.OK;

  }

  public @Nonnull MZmineModuleCategory getModuleCategory() {
    return MZmineModuleCategory.DATAANALYSIS;
  }

  @Override
  public @Nonnull Class<? extends ParameterSet> getParameterSetClass() {
    return MummichogParameters.class;
  }
}
