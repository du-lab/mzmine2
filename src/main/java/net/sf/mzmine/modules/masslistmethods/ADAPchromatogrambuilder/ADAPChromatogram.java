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
 *
 * Edited and modified by Owen Myers (Oweenm@gmail.com)
 */


package net.sf.mzmine.modules.masslistmethods.ADAPchromatogrambuilder;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.SortedMap;
import java.util.TreeMap;
import javax.annotation.Nonnull;
import com.google.common.collect.Range;
import com.google.common.primitives.Ints;
import net.sf.mzmine.datamodel.DataPoint;
import net.sf.mzmine.datamodel.Feature;
import net.sf.mzmine.datamodel.IsotopePattern;
import net.sf.mzmine.datamodel.RawDataFile;
import net.sf.mzmine.datamodel.Scan;
import net.sf.mzmine.datamodel.impl.SimplePeakInformation;
import net.sf.mzmine.main.MZmineCore;
import net.sf.mzmine.util.ScanUtils;



/**
 * Chromatogram implementing ChromatographicPeak.
 */
public class ADAPChromatogram implements Feature {
  private SimplePeakInformation peakInfo;

  // Data file of this chromatogram
  private RawDataFile dataFile;

  private int[] scanNumbers;
  private DataPoint[] dataPoints;

  // Data points of the chromatogram (map of scan number -> m/z peak)
  // This variable will be destroyed in finishChromatogram()
  private SortedMap<Integer, DataPoint> dataPointsMap;

  // Chromatogram m/z, RT, height, area. The mz value will be the highest points mz value
  private double mz, rt, height, area, weightedMz;
  private Double fwhm = null, tf = null, af = null;

  // Top intensity scan, fragment scan
  private int representativeScan = -1, fragmentScan = -1;

  // All MS2 fragment scan numbers
  private int[] allMS2FragmentScanNumbers = new int[] {};

  // Ranges of raw data points
  private Range<Double> rawDataPointsIntensityRange, rawDataPointsMZRange, rawDataPointsRTRange;

  // Isotope pattern. Null by default but can be set later by deisotoping
  // method.
  private IsotopePattern isotopePattern;
  private int charge = 0;

  // Victor Trevino
  private double mzSum = 0;
  private int mzN = 0;
  private double weightedMzSum = 0;
  private int weightedMzN = 0;
  private double sumOfWeights = 0;

  // Remove these variables when chromatogram is finished
  final private double minMz;
  final private double maxMz;

  /**
   * Initializes this Chromatogram
   */
  public ADAPChromatogram(RawDataFile dataFile, double minMz, double maxMz) {

    this.dataFile = dataFile;
    this.minMz = minMz;
    this.maxMz = maxMz;
    rawDataPointsRTRange = dataFile.getDataRTRange(1);

    dataPointsMap = new TreeMap<>();
  }

  public int findNumberOfContinuousPointsAboveNoise(double noise, int[] scanNumbers) {
    // sort the array containing all of the scan numbers of the point added
    // loop over the sorted array now.
    // if you find a point with intensity higher than noise start the count
    // if the next scan contains a point higher than the noise update the count
    // otherwise start it oer when you hit a sufficiently high point.
    // keep track of the largest count which will be returned.

    int[] sortedScanNumbers = Arrays.copyOf(scanNumbers, scanNumbers.length);
    Arrays.sort(sortedScanNumbers);

    int bestCount = 0;
    int curCount = 0;
    for (int scanNumber : sortedScanNumbers) {

      DataPoint dataPoint = getDataPoint(scanNumber);

      if (dataPoint == null) {
        curCount = 0;
        continue;
      }

      if (dataPoint.getIntensity() <= noise) {
        curCount = 0;
        continue;
      }

      ++curCount;

      if (curCount > bestCount)
        bestCount = curCount;
    }
    return bestCount;
  }

  /**
   * This method adds a MzPeak to this Chromatogram. All values of this Chromatogram (rt, m/z,
   * intensity and ranges) are updated on request
   * 
   * @param mzValue
   */
  public void addMzPeak(int scanNumber, DataPoint mzValue) {

    if (mzValue.getMZ() < minMz || mzValue.getMZ() >= maxMz)
      throw new IllegalStateException("Attempt to add m/z-value outside chromatogram's m/z-range");

    double curIntensity;

    // If we already have a mzvalue for the scan number then we need to add the intensities
    // together before putting it into the dataPointsMap, otherwise the chromatogram is only
    // representing the intesities of the llast added point for that scan.
    //
    // For now just don't add the point if we have it already. The highest point will be the
    // first one added
    if (dataPointsMap.containsKey(scanNumber))
      return;

    dataPointsMap.put(scanNumber, mzValue);
    mzSum += mzValue.getMZ();
    mzN++;
    mz = mzSum / mzN;

    curIntensity = mzValue.getIntensity();

    weightedMzN++;
    weightedMzSum += curIntensity * mzValue.getMZ();
    sumOfWeights += curIntensity;

    weightedMz = weightedMzSum / sumOfWeights;
  }

  @Override
  public DataPoint getDataPoint(int scanNumber) throws IllegalStateException {

    if (dataPointsMap != null)
      return dataPointsMap.get(scanNumber);

    if (scanNumbers != null) {
      int index = Arrays.binarySearch(scanNumbers, scanNumber);
      return (index < 0) ? null : dataPoints[index];
    }

    throw new IllegalStateException("Both dataPointsMap and scanNumbers are null.");
  }

  /**
   * This method returns m/z value of the chromatogram
   */
  @Override
  public double getMZ() {
    return mz;
  }

  /**
   * This method returns weighted mean of m/z values comprising the chromatogram
   */
  public double getWeightedMZ() {
    return weightedMz;
  }


  /**
   * This method returns a string with the basic information that defines this peak
   * 
   * @return String information
   */
  @Override
  public String toString() {
    return "Chromatogram " + MZmineCore.getConfiguration().getMZFormat().format(mz) + " m/z";
  }

  @Override
  public double getArea() {
    return area;
  }

  @Override
  public double getHeight() {
    return height;
  }

  @Override
  public int getMostIntenseFragmentScanNumber() {
    return fragmentScan;
  }

  @Override
  public int[] getAllMS2FragmentScanNumbers() {
    return allMS2FragmentScanNumbers;
  }

  @Override
  public @Nonnull FeatureStatus getFeatureStatus() {
    return FeatureStatus.DETECTED;
  }

  @Override
  public double getRT() {
    return rt;
  }

  @Override
  public @Nonnull Range<Double> getRawDataPointsIntensityRange() {
    return rawDataPointsIntensityRange;
  }

  @Override
  public @Nonnull Range<Double> getRawDataPointsMZRange() {
    return rawDataPointsMZRange;
  }

  @Override
  public @Nonnull Range<Double> getRawDataPointsRTRange() {
    return rawDataPointsRTRange;
  }

  @Override
  public int getRepresentativeScanNumber() {
    return representativeScan;
  }

  @Override
  public @Nonnull int[] getScanNumbers() {
    return scanNumbers;
  }

  @Override
  public @Nonnull RawDataFile getDataFile() {
    return dataFile;
  }

  @Override
  public IsotopePattern getIsotopePattern() {
    return isotopePattern;
  }

  @Override
  public void setIsotopePattern(@Nonnull IsotopePattern isotopePattern) {
    this.isotopePattern = isotopePattern;
  }

  @Override
  public void outputChromToFile() {

    int[] allScanNumbers = getScanNumbers();
    Arrays.sort(allScanNumbers);
    try {
      String fileName = dataFile.getName();

      String fileNameNoExtention = fileName.split("\\.")[0];

      // PrintWriter fileWriter = new PrintWriter("owen_mzmine_chrom_out.txt","UTF-8");
      PrintWriter fileWriter = new PrintWriter(
          new FileOutputStream(new File(fileNameNoExtention + "_mzmine_chrom_out.csv"), true));

      String mzStr = String.valueOf(getMZ());
      // fileWriter.println("New chromatogram. mz: "+mzStr);
      String curIntStr;
      String curRtStr;
      String curScStr;
      String curMzStr;
      double curInt;
      double curRt;
      double curMz;
      String spacer = ",";


      for (int i : allScanNumbers) {


        curRt = dataFile.getScan(i).getRetentionTime();
        DataPoint mzPeak = getDataPoint(i);

        if (mzPeak == null) {
          curInt = 0.0;
          curMz = -1;
        } else {
          curInt = mzPeak.getIntensity();
          curMz = mzPeak.getMZ();
        }

        curIntStr = String.valueOf(curInt);
        curRtStr = String.valueOf(curRt);
        curScStr = String.valueOf(i);
        curMzStr = String.valueOf(curMz);

        fileWriter.println(curRtStr + spacer + curIntStr + spacer + curMzStr);
      }
      fileWriter.close();
    }

    catch (FileNotFoundException ex) {
      ex.printStackTrace();
    }

  }


  public void finishChromatogram() {

    // Update raw data point ranges, height, rt and representative scan
    height = Double.MIN_VALUE;

    double minMz = mz;
    double maxMz = mz;
    double minRetTime = Double.MAX_VALUE;
    double maxRetTime = -Double.MAX_VALUE;
    double minIntensity = Double.MAX_VALUE;
    double maxIntensity = -Double.MAX_VALUE;
    for (Map.Entry<Integer, DataPoint> e : dataPointsMap.entrySet()) {

      DataPoint mzPeak = e.getValue();
      int currentScan = e.getKey();
      double currentIntensity = mzPeak.getIntensity();
      double currentMz = mzPeak.getMZ();
      double currentRetTime = dataFile.getScan(e.getKey()).getRetentionTime();

      if (height < currentIntensity) {
        height = currentIntensity;
        rt = currentRetTime;
        mz = currentMz;
        representativeScan = currentScan;
      }

      if (currentMz < minMz) minMz = currentMz;
      if (currentMz > maxMz) maxMz = currentMz;

      if (currentIntensity < minIntensity) minIntensity = currentIntensity;
      if (currentIntensity > maxIntensity) maxIntensity = currentIntensity;

      if (currentRetTime > maxRetTime) maxRetTime = currentRetTime;
      if (currentRetTime < minRetTime) minRetTime = currentRetTime;
    }

    rawDataPointsMZRange = Range.closed(minMz, maxMz);
    rawDataPointsIntensityRange = Range.closed(minIntensity, maxIntensity);
    rawDataPointsRTRange = Range.closed(minRetTime, maxRetTime);

    int[] allScanNumbers = Ints.toArray(dataPointsMap.keySet());
    Arrays.sort(allScanNumbers);

    // Update area
    area = 0;
    for (int i = 1; i < allScanNumbers.length; i++) {
      // For area calculation, we use retention time in seconds
      double previousRT = dataFile.getScan(allScanNumbers[i - 1]).getRetentionTime() * 60d;
      double currentRT = dataFile.getScan(allScanNumbers[i]).getRetentionTime() * 60d;
      double previousHeight = dataPointsMap.get(allScanNumbers[i - 1]).getIntensity();
      double currentHeight = dataPointsMap.get(allScanNumbers[i]).getIntensity();
      area += (currentRT - previousRT) * (currentHeight + previousHeight) / 2;
    }

    // Update fragment scan
    fragmentScan =
        ScanUtils.findBestFragmentScan(dataFile, dataFile.getDataRTRange(1), rawDataPointsMZRange);

    allMS2FragmentScanNumbers = ScanUtils.findAllMS2FragmentScans(dataFile,
        dataFile.getDataRTRange(1), rawDataPointsMZRange);

    if (fragmentScan > 0) {
      Scan fragmentScanObject = dataFile.getScan(fragmentScan);
      int precursorCharge = fragmentScanObject.getPrecursorCharge();
      if (precursorCharge > 0)
        this.charge = precursorCharge;
    }

    // Save DataPoints from dataPointsMap to scanNumbers[] and dataPoints[]
    this.scanNumbers = new int[dataPointsMap.size()];
    this.dataPoints = new DataPoint[dataPointsMap.size()];
    int i = 0;
    for (SortedMap.Entry<Integer, DataPoint> e : dataPointsMap.entrySet()) {
      this.scanNumbers[i] = e.getKey();
      this.dataPoints[i] = e.getValue();
      ++i;
    }

    // Discard the fields we don't need anymore
    dataPointsMap.clear();
    dataPointsMap = null;
  }

  @Override
  public int getCharge() {
    return charge;
  }

  @Override
  public void setCharge(int charge) {
    this.charge = charge;
  }

  @Override
  public Double getFWHM() {
    return fwhm;
  }

  @Override
  public void setFWHM(Double fwhm) {
    this.fwhm = fwhm;
  }

  @Override
  public Double getTailingFactor() {
    return tf;
  }

  @Override
  public void setTailingFactor(Double tf) {
    this.tf = tf;
  }

  @Override
  public Double getAsymmetryFactor() {
    return af;
  }

  @Override
  public void setAsymmetryFactor(Double af) {
    this.af = af;
  }

  @Override
  public void setPeakInformation(SimplePeakInformation peakInfoIn) {
    this.peakInfo = peakInfoIn;
  }

  @Override
  public SimplePeakInformation getPeakInformation() {
    return peakInfo;
  }

  public boolean inMzRange(double mz) {
    return minMz <= mz && mz < maxMz;
  }

  public Double getMaxMz() {
    return maxMz;
  }

  public Double getMinMz() {
    return minMz;
  }

}
