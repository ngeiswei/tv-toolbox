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

-- Plot histogram
plotMaps :: Bool -> Bool -> String -> [(String, (Map Double Double))] -> IO ()
plotMaps save zoom title names_maps =
  plotPathsStyle attributes (Prelude.map fmt names_maps)
  where attributes = [Title title, XLabel "Probability", YLabel "Density"]
                     ++ (if zoom then [] else [XRange (0.0, 1.0)])
                     ++ (if save then [PNG (title ++ ".png")] else [])
        fmt (n, m) = (defaultStyle {lineSpec = CustomStyle [LineTitle n]},
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
    alpha = 2
    beta = 5
    smp_size = 1000000 -- Number of samples
    nbr_bins = 3       -- NEXT: when nbr_bins is low what to do?
                       -- Should the bin function be shifted by
                       -- half-bin to the right?
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
    False False
    (format "Sample vs Fit (alpha={0}, beta={1}, bins={2})"
     [show alpha, show beta, show nbr_bins])
    [((format "Sample (size={0})" [show smp_size]), smp_map),
     ((format "Fitted (alpha={0}, beta={1})" [show smp_alpha, show smp_beta]),
      fit_map)]
  print (format
         ("Beta Distribution alpha = {0}, beta = {1}, cdf_half = {2}"
          ++ ", normalized_histo = {3}"
          ++ ", smp_mean = {4}, smp_var = {5}, smp_stdev = {6}"
          ++ ", smp_alpha = {7}, smp_beta = {8}")
         [show (bdAlpha bd), show (bdBeta bd), show (cdf_half),
          show (histogramToPdfMap nbr_bins smp_histo),
          show smp_mean, show smp_var, show smp_stdev,
          show smp_alpha, show smp_beta])
