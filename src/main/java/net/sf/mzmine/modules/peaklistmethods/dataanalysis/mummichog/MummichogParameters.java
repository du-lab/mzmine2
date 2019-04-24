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

import java.text.DecimalFormat;
import net.sf.mzmine.parameters.Parameter;
import net.sf.mzmine.parameters.impl.SimpleParameterSet;
import net.sf.mzmine.parameters.parametertypes.BooleanParameter;
import net.sf.mzmine.parameters.parametertypes.ComboParameter;
import net.sf.mzmine.parameters.parametertypes.DoubleParameter;
import net.sf.mzmine.parameters.parametertypes.StringParameter;
import net.sf.mzmine.parameters.parametertypes.filenames.DirectoryParameter;
import net.sf.mzmine.parameters.parametertypes.selectors.PeakListsParameter;

public class MummichogParameters extends SimpleParameterSet {

  public static final PeakListsParameter peakLists = new PeakListsParameter(1, 1);

  public static final DoubleParameter cutoff = new DoubleParameter("Significance Cutoff",
      "Significance cutoff p-value in input file", DecimalFormat.getInstance(), 0.05);
  public static final ComboParameter<String> network =
      new ComboParameter<String>("Metabolic Network", "Choose the metabolic network to be used",
          new String[] {"human", "worm"}, "human");
  public static final ComboParameter<String> modeling = new ComboParameter<String>(
      "Distriution Estimator", "Choose the type of estimation to be used",
      new String[] {"Non-Parametric", "Gamma"}, "Non-Parametric");

  public static final BooleanParameter force_primary_ion = new BooleanParameter(
      "Primary Ion Mandate", "Enforce Primary Ion in Empirical Compunds", true);

  public static final DirectoryParameter output =
      new DirectoryParameter("Output Directory", "Output directory");

  public MummichogParameters() {
    super(new Parameter[] {peakLists, cutoff, network, force_primary_ion, modeling, output});
  }
}
