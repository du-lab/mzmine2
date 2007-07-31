/*
 * Copyright 2006-2007 The MZmine Development Team
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

package net.sf.mzmine.modules.dataanalysis.intensityplot;

import java.text.Format;

import net.sf.mzmine.data.Peak;
import net.sf.mzmine.io.RawDataFile;
import net.sf.mzmine.main.MZmineCore;
import net.sf.mzmine.userinterface.Desktop;

import org.jfree.chart.labels.CategoryToolTipGenerator;
import org.jfree.chart.labels.XYToolTipGenerator;
import org.jfree.data.category.CategoryDataset;
import org.jfree.data.xy.XYDataset;

/**
 * 
 */
class IntensityPlotTooltipGenerator implements CategoryToolTipGenerator,
        XYToolTipGenerator {

    /**
     * @see org.jfree.chart.labels.CategoryToolTipGenerator#generateToolTip(org.jfree.data.category.CategoryDataset,
     *      int, int)
     */
    public String generateToolTip(CategoryDataset dataset, int row, int column) {
        Desktop desktop = MZmineCore.getDesktop();
        Format intensityFormat = desktop.getIntensityFormat();
        Peak peaks[] = ((IntensityPlotDataset) dataset).getPeaks(row, column);
        RawDataFile files[] = ((IntensityPlotDataset) dataset).getFiles(column);

        StringBuffer sb = new StringBuffer();
        sb.append("<html>");
        for (int i = 0; i < files.length; i++) {
            sb.append(files[i].toString());
            sb.append(": ");
            if (peaks[i] != null) {
                sb.append(peaks[i].toString());
                sb.append(", height: ");
                sb.append(intensityFormat.format(peaks[i].getHeight()));
            } else {
                sb.append("N/A");
            }
            sb.append("<br>");
        }

        sb.append("</html>");

        return sb.toString();
    }

    public String generateToolTip(XYDataset dataset, int series, int item) {
        return generateToolTip((CategoryDataset) dataset, series, item);
    }

}
