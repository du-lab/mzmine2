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

package net.sf.mzmine.util;

import java.text.Format;

import net.sf.mzmine.data.CompoundIdentity;
import net.sf.mzmine.data.IsotopePattern;
import net.sf.mzmine.data.Peak;
import net.sf.mzmine.data.PeakListRow;
import net.sf.mzmine.main.MZmineCore;

/**
 * Utilities for peaks and peak lists
 * 
 */
public class PeakUtils {

    /**
     * Common utility method to be used as Peak.toString() method in various
     * Peak implementations
     * 
     * @param peak Peak to be converted to String
     * @return String representation of the peak
     */
    public static String peakToString(Peak peak) {
        StringBuffer buf = new StringBuffer();
        Format mzFormat = MZmineCore.getMZFormat();
        Format timeFormat = MZmineCore.getRTFormat();
        buf.append(mzFormat.format(peak.getMZ()));
        buf.append(" m/z @");
        buf.append(timeFormat.format(peak.getRT()));
        return buf.toString();
    }

    /**
     * Returns peak with the lowest m/z value of the isotope pattern
     */
    public static Peak getLowestMZPeak(IsotopePattern pattern) {

        if (pattern == null)
            return null;

        Peak[] peaks = pattern.getOriginalPeaks();

        if ((peaks == null) || (peaks.length == 0))
            return null;

        Peak lowestMZPeak = peaks[0];
        for (Peak peak : peaks)
            if (peak.getMZ() < lowestMZPeak.getMZ())
                lowestMZPeak = peak;

        return lowestMZPeak;

    }

    /**
     * Returns the most intense peak of the isotope pattern
     */
    public static Peak getMostIntensePeak(IsotopePattern pattern) {

        if (pattern == null)
            return null;

        Peak[] peaks = pattern.getOriginalPeaks();

        if ((peaks == null) || (peaks.length == 0))
            return null;

        Peak mostIntensePeak = peaks[0];
        for (Peak peak : peaks)
            if (peak.getArea() > mostIntensePeak.getArea())
                mostIntensePeak = peak;

        return mostIntensePeak;

    }

    /**
     * Returns average m/z of peaks and/or isotope patterns on the peak list
     * row. For isotope patterns, uses the lowest m/z value of the pattern.
     */
    public static float getAverageMZUsingLowestMZPeaks(PeakListRow peakListRow) {

        if (peakListRow == null)
            return 0.0f;

        Peak[] peaks = peakListRow.getPeaks();

        if ((peaks == null) || (peaks.length == 0))
            return 0.0f;

        float mzSum = 0.0f;
        for (Peak peak : peaks) {
            if (peak instanceof IsotopePattern)
                mzSum += getLowestMZPeak((IsotopePattern) peak).getMZ();
            else
                mzSum += peak.getMZ();
        }

        return mzSum / (float) peaks.length;

    }

    /**
     * Returns average m/z of the isotope patterns and peaks on the row. For
     * isotope patterns, uses m/z value of the most intense peak of the pattern.
     */
    public static float getAverageMZUsingMostIntensePeaks(
            PeakListRow peakListRow) {

        if (peakListRow == null)
            return 0.0f;

        Peak[] peaks = peakListRow.getPeaks();

        if ((peaks == null) || (peaks.length == 0))
            return 0.0f;

        float mzSum = 0.0f;
        for (Peak peak : peaks) {
            if (peak instanceof IsotopePattern)
                mzSum += getMostIntensePeak((IsotopePattern) peak).getMZ();
            else
                mzSum += peak.getMZ();
        }

        return mzSum / (float) peaks.length;

    }

    /**
     * Compares identities of two peak list rows. 1) if preferred identities are
     * available, they must be same 2) if no identities are available on both
     * rows, return true 3) otherwise all identities on both rows must be same
     * 
     * @return True if identities match between rows
     * 
     * TODO Compare using IDs instead of names when possible
     */
    public static boolean compareIdentities(PeakListRow row1, PeakListRow row2) {

        if ((row1 == null) || (row2 == null))
            return false;

        // If both have preferred identity available, then compare only those
        CompoundIdentity row1PreferredIdentity = row1.getPreferredCompoundIdentity();
        CompoundIdentity row2PreferredIdentity = row2.getPreferredCompoundIdentity();
        if ((row1PreferredIdentity != null) && (row2PreferredIdentity != null)) {
            if (row1PreferredIdentity.getCompoundName().equals(
                    row2PreferredIdentity.getCompoundName()))
                return true;
            else
                return false;
        }

        // If no identities at all for both rows, then return true
        CompoundIdentity[] row1Identities = row1.getCompoundIdentities();
        CompoundIdentity[] row2Identities = row2.getCompoundIdentities();
        if ((row1Identities.length == 0) && (row2Identities.length == 0))
            return true;

        // Otherwise compare all against all and require that each identity has
        // a matching identity on the other row
        if (row1Identities.length != row2Identities.length)
            return false;
        boolean sameID = false;
        for (CompoundIdentity row1Identity : row1Identities) {
            sameID = false;
            for (CompoundIdentity row2Identity : row2Identities) {
                if (row1Identity.getCompoundName().equals(
                        row2Identity.getCompoundName())) {
                    sameID = true;
                    break;
                }
            }
            if (!sameID)
                break;
        }

        return sameID;

    }

    /**
     * Returns true if peak list row contains a compound identity matching to id
     * 
     * TODO Compare IDs instead of names when possible
     */
    public static boolean containsIdentity(PeakListRow row, CompoundIdentity id) {

        for (CompoundIdentity identity : row.getCompoundIdentities()) {
            if (identity.getCompoundName().equals(id.getCompoundName()))
                return true;
        }

        return false;
    }

}
