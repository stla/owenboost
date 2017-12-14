{-# LANGUAGE CPP, ForeignFunctionInterface #-}
{-# LANGUAGE QuasiQuotes              #-}
{-# LANGUAGE TemplateHaskell          #-}

module Owen
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

studentC :: CDouble -> CInt -> [CDouble] -> IO [CDouble]
studentC q nu delta = do
  mdelta <- V.thaw (V.fromList delta)
  ptr <- [CPP.exp| double* {
      studentC($(double q), $(int nu), $vec-ptr:(double* mdelta), $vec-len:mdelta)
    } |]
  VM.clear mdelta
  peekArray (length delta) ptr

studentCDF :: CDouble -> CInt -> [CDouble] -> IO [CDouble]
studentCDF q nu delta = do
  mdelta <- V.thaw (V.fromList delta)
  ptr <- [CPP.exp| double* {
      studentCDF($(double q), $(int nu), $vec-ptr:(double* mdelta), $vec-len:mdelta)
    } |]
  peekArray (length delta) ptr

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
