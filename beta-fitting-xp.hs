#!/usr/bin/env runhaskell

-- Attempt to fit a distribution with a beta distribution

import Control.Monad (replicateM)
import System.Random.MWC (initialize)
import Data.Map (Map, map, toList, fromList)
import Data.Vector.Unboxed (singleton)
import Data.Histogram (Histogram, fromList, toList, toMap, size)
import Statistics.Distribution.Beta (BetaDistribution, betaDistr, bdAlpha, bdBeta)
import Statistics.Distribution (cumulative, density, genContVar)
import Text.Format (format)
import Graphics.Gnuplot.Simple (plotPathStyle, plotPathsStyle,
                                Attribute(Title, XLabel, YLabel, XRange, PNG),
                                PlotStyle, defaultStyle,
                                plotType, PlotType(Lines, Steps, Boxes),
                                lineSpec, LineSpec(CustomStyle),
                                LineAttr(LineTitle))

import TVToolBox (bin, square)

---------------
-- Functions --
---------------

-- Take an histogram and output a Map Double Double where the counts
-- have been normalized to sum up to 1.
histogramToPdfMap :: Integer -> (Histogram Double) -> (Map Double Double)
histogramToPdfMap nbr_bins hist =
  (Data.Map.map (\v -> ((fromIntegral v) * c)) (toMap hist))
  where s = fromIntegral (size hist)
        n = fromIntegral nbr_bins
        c = n / s

-- Plot histogram.  Take a list of GnuPlot attributes, a title and a
-- list of triple (name, plot type, map).  It automatically add a
-- "Probability" XLabel and a "Density" YLabel to the list of
-- attributes.
plotMaps :: [Attribute] -> String -> [(String, PlotType, (Map Double Double))] -> IO ()
plotMaps attributes title names_types_maps =
  plotPathsStyle attributes_xt (Prelude.map fmt names_types_maps)
  where attributes_xt = [Title title, XLabel "Probability", YLabel "Density"]
                        ++ attributes
        fmt (n, t, m) = (defaultStyle {plotType = t,
                                       lineSpec = CustomStyle [LineTitle n]},
                         (Data.Map.toList m))

-- Produce histogram with n bins from a sample
sampleToHistogram :: Integer -> [Double] -> (Histogram Double)
sampleToHistogram n smp = Data.Histogram.fromList (Prelude.map (bin n) smp)

-- Estimate mean from sample
sampleToMean :: [Double] -> Double
sampleToMean smp = (sum smp) / (fromIntegral (length smp))

-- Estimate variance from sample
sampleToVariance :: [Double] -> Double
sampleToVariance smp = sampleToMean (Prelude.map (\x -> square (x - smp_mean)) smp)
  where smp_mean = sampleToMean smp

-- Estimate alpha parameter of beta distribution from sample
sampleToAlpha :: [Double] -> Double
sampleToAlpha smp = (m*(1 - m)/v - 1)*m
  where m = sampleToMean smp
        v = sampleToVariance smp

-- Estimate beta parameter of beta distribution from sample
sampleToBeta :: [Double] -> Double
sampleToBeta smp = (m*(1 - m)/v - 1)*(1 - m)
  where m = sampleToMean smp
        v = sampleToVariance smp

-- Produce a map representing the PDF of a beta distribution
-- discretizing with n bins.
bdToPdfMap :: Integer -> BetaDistribution -> (Map Double Double)
bdToPdfMap n bd = (Data.Map.fromList [(toProb i, toDensity i) | i <- [0..n]])
  where
    nd = (fromIntegral n) :: Double
    toProb i = (fromIntegral i) / nd
    toDensity i = density bd (toProb i)

----------
-- Main --
----------

main :: IO ()
main = do
  let
    -- Constants
    seed = 0
    alpha = 20
    beta = 5
    smp_size = 10000 -- Number of samples
    nbr_bins = 40
    -- Define beta distribution
    bd = betaDistr alpha beta
    cdf_half = cumulative bd 0.5
  -- Initialize random generator
  g <- initialize (singleton seed)
  -- Sample beta distribution
  smp <- replicateM smp_size (genContVar bd g)
  let
    -- Produce histogram from sample
    smp_histo = sampleToHistogram nbr_bins smp
    -- Produce normalized map for plotting
    smp_map = histogramToPdfMap nbr_bins smp_histo
    -- Find corresponding alpha and beta
    smp_mean = sampleToMean smp
    smp_var = sampleToVariance smp
    smp_stdev = sqrt smp_var
    smp_alpha = sampleToAlpha smp
    smp_beta = sampleToBeta smp
    -- Fit beta distribution
    fit_bd = betaDistr smp_alpha smp_beta
    -- Produce PDF map for plotting
    fit_map = bdToPdfMap nbr_bins fit_bd
  plotMaps
    []
    (format "Sampled vs Fitted (alpha={0}, beta={1}, bins={2})"
     [show alpha, show beta, show nbr_bins])
    [((format "Sampled (size={0})" [show smp_size]), Boxes, smp_map),
     ((format "Fitted (alpha={0}, beta={1})" [show smp_alpha, show smp_beta]),
      Lines, fit_map)]
  print (format
         ("Beta Distribution alpha = {0}, beta = {1}, cdf_half = {2}"
          ++ ", normalized_histo = {3}"
          ++ ", smp_mean = {4}, smp_var = {5}, smp_stdev = {6}"
          ++ ", smp_alpha = {7}, smp_beta = {8}")
         [show (bdAlpha bd), show (bdBeta bd), show (cdf_half),
          show (histogramToPdfMap nbr_bins smp_histo),
          show smp_mean, show smp_var, show smp_stdev,
          show smp_alpha, show smp_beta])
