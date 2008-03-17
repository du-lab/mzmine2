/*
 * Copyright 2006-2008 The MZmine Development Team
 * 
 * This file is part of MZmine.
 * 
 * MZmine is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 * 
 * MZmine is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * MZmine; if not, write to the Free Software Foundation, Inc., 51 Franklin St,
 * Fifth Floor, Boston, MA 02110-1301 USA
 */

package net.sf.mzmine.modules.visualization.twod;

import net.sf.mzmine.data.Parameter;
import net.sf.mzmine.data.ParameterType;
import net.sf.mzmine.data.impl.SimpleParameter;
import net.sf.mzmine.data.impl.SimpleParameterSet;
import net.sf.mzmine.main.MZmineCore;

/**
 * 2D visualizer parameter set
 */
public class TwoDParameters extends SimpleParameterSet {

    public static final Integer[] msLevels = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };

    public static final Parameter msLevel = new SimpleParameter(
            ParameterType.INTEGER, "MS level", "MS level of plotted scans", 1,
            msLevels);

    public static final Parameter minRT = new SimpleParameter(
            ParameterType.FLOAT, "Minimum retention time",
            "X axis lower limit", null, new Float(0.0), new Float(0.0), null,
            MZmineCore.getRTFormat());

    public static final Parameter maxRT = new SimpleParameter(
            ParameterType.FLOAT, "Maximum retention time",
            "X axis upper limit", null, new Float(600.0), new Float(0.0), null,
            MZmineCore.getRTFormat());

    public static final Parameter minMZ = new SimpleParameter(
            ParameterType.FLOAT, "Minimum M/Z", "m/z lower limit", "m/z",
            new Float(100.0), new Float(0.0), null, MZmineCore.getMZFormat());

    public static final Parameter maxMZ = new SimpleParameter(
            ParameterType.FLOAT, "Maximum M/Z", "m/z upper limit", "m/z",
            new Float(1000.0), new Float(0.0), null, MZmineCore.getMZFormat());

    public TwoDParameters() {
        super(new Parameter[] { msLevel, minRT, maxRT, minMZ, maxMZ });
    }

}
