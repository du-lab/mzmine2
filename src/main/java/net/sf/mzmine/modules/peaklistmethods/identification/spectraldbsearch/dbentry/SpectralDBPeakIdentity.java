/*
 * Copyright 2006-2018 The MZmine 2 Development Team
 * 
 * This file is part of MZmine 2.
 * 
 * MZmine 2 is free software; you can redistribute it and/or modify it under the terms of the GNU
 * General Public License as published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 * 
 * MZmine 2 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
 * even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with MZmine 2; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
 * USA
 */

package net.sf.mzmine.modules.peaklistmethods.identification.spectraldbsearch.dbentry;

import java.text.DecimalFormat;
import java.text.MessageFormat;
import net.sf.mzmine.datamodel.impl.SimplePeakIdentity;
import net.sf.mzmine.util.maths.similarity.SpectraSimilarity;

public class SpectralDBPeakIdentity extends SimplePeakIdentity {
  private static final DecimalFormat COS_FORM = new DecimalFormat("0.000");

  private final SpectralDBEntry entry;
  private final SpectraSimilarity similarity;

  public SpectralDBPeakIdentity(SpectralDBEntry entry, SpectraSimilarity similarity,
      String method) {
    super(
        MessageFormat.format("{0} as {3} ({1}) {2} cos={4}",
            entry.getField(DBEntryField.NAME).orElse("NONAME"), entry.getPrecursorMZ(),
            entry.getField(DBEntryField.FORMULA).orElse(""),
            entry.getField(DBEntryField.IONTYPE).orElse("N/A"),
            COS_FORM.format(similarity.getCosine())),
        entry.getField(DBEntryField.FORMULA).orElse("").toString(), method, "", "");
    this.entry = entry;
    this.similarity = similarity;
  }

  public SpectralDBEntry getEntry() {
    return entry;
  }

  public SpectraSimilarity getSimilarity() {
    return similarity;
  }

}
