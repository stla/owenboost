{-# LANGUAGE CPP, ForeignFunctionInterface #-}
{-# LANGUAGE QuasiQuotes              #-}
{-# LANGUAGE TemplateHaskell          #-}

module Test
  where
import           Data.Monoid                  ((<>))
import qualified Data.Vector.Storable         as V
import qualified Data.Vector.Storable.Mutable as VM
import           Foreign
import           Foreign.C.Types
import qualified Language.C.Inline            as C
import qualified Language.C.Inline.Cpp        as CPP

CPP.context (CPP.cppCtx <> C.vecCtx)

CPP.include "<owen128.cpp>"

sumVec :: VM.IOVector CDouble -> IO CDouble
sumVec vec = [CPP.block| double {
    double sum = 0;
    int i;
    for (i = 0; i < $vec-len:vec; i++) {
      sum += $vec-ptr:(double *vec)[i];
    }
    return sum;
  } |]

testvec :: IO ()
testvec = do
  x <- sumVec =<< V.thaw (V.fromList [1,2,3])
  print x

-- doubleVec :: [CDouble] -> IO [CDouble]
-- doubleVec vec = do
--   mvec <- V.thaw (V.fromList vec)
--   ptr <- [CPP.block| double* {
--     double *out = new double[$vec-len:mvec];
--     int i;
--     for (i = 0; i < $vec-len:mvec; i++) {
--       out[i] = 2 * $vec-ptr:(double *mvec)[i];
--     }
--     return out;
--   } |]
--   peekArray (length vec) ptr

doubleVec2 :: [CDouble] -> IO [CDouble]
doubleVec2 x = do
  xx <- V.thaw (V.fromList x)
  ptr <- [CPP.exp| double* { testvectorout($vec-ptr:(double* xx), $vec-len:xx) } |]
  peekArray (length x) ptr

studentC :: CDouble -> CInt -> [CDouble] -> IO [CDouble]
studentC q nu delta = do
  mdelta <- V.thaw (V.fromList delta)
  ptr <- [CPP.exp| double* {
      studentC($(double q), $(int nu), $vec-ptr:(double* mdelta), $vec-len:mdelta)
    } |]
  VM.clear mdelta
  peekArray (length delta) ptr

pnorm :: CDouble
pnorm = [CPP.pure| double { pnorm(2.0) }|]

studentCDF :: CDouble -> CInt -> [CDouble] -> IO [CDouble]
studentCDF q nu delta = do
  mdelta <- V.thaw (V.fromList delta)
  ptr <- [CPP.exp| double* {
      studentCDF($(double q), $(int nu), $vec-ptr:(double* mdelta), $vec-len:mdelta)
    } |]
  peekArray (length delta) ptr

f :: Double
f = 2

owenC :: CInt -> CDouble -> [CDouble] -> [CDouble] -> IO [CDouble]
owenC nu t delta r = do
  mdelta <- V.thaw (V.fromList delta)
  mr <- V.thaw (V.fromList r)
  ptr <- [CPP.exp| double* {
      owenC($(int nu), $(double t), $vec-ptr:(double* mdelta), $vec-ptr:(double* mr), $vec-len:mdelta)
    } |]
  VM.clear mdelta
  VM.clear mr
  peekArray (length delta) ptr

owenQ :: CInt -> CDouble -> [CDouble] -> [CDouble] -> IO [CDouble]
owenQ nu t delta r = do
  mdelta <- V.thaw (V.fromList delta)
  mr <- V.thaw (V.fromList r)
  ptr <- [CPP.exp| double* {
      owenQ($(int nu), $(double t), $vec-ptr:(double* mdelta), $vec-ptr:(double* mr), $vec-len:mdelta)
    } |]
  VM.clear mdelta
  VM.clear mr
  peekArray (length delta) ptr

owenQ128 :: CInt -> CDouble -> [CDouble] -> [CDouble] -> IO [CDouble]
owenQ128 nu t delta r = do
  mdelta <- V.thaw (V.fromList delta)
  mr <- V.thaw (V.fromList r)
  ptr <- [CPP.exp| double* {
      owenQ128($(int nu), $(double t), $vec-ptr:(double* mdelta), $vec-ptr:(double* mr), $vec-len:mdelta)
    } |]
  VM.clear mdelta
  VM.clear mr
  peekArray (length delta) ptr
